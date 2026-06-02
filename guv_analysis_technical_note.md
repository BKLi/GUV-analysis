# GUV Fluorescence Analysis Pipeline — Technical Note

**Script:** `guv_analysis.ipynb`  
**Date:** June 2026  
**Author:** Bingkun Li

---

## 1 · Overview and pipeline logic

The pipeline segments Giant Unilamellar Vesicles (GUVs) from two-channel confocal fluorescence images (GFP + mCherry) stored as CellSens `.vsi` files, then measures fluorescence intensity independently in the **membrane ring** and **lumen** of each GUV. The core analytical chain is:

```
VSI file
  └─ read_vsi()            → (C, Y, X) uint16 array
       └─ split_channels() → gfp_img, mcherry_img  (2-D uint16)
            └─ normalise() → 0–255 float32  (for Cellpose input)
                 └─ segment_guvs()          → integer label mask  (int32)
                      └─ [napari correction] → masks_corrected
                           └─ make_membrane_lumen_masks() → membrane/lumen dicts
                                └─ quantify()  (×2, one per channel)
                                     └─ DataFrame  (per GUV rows)
                                          └─ _bulk_mode()  → background subtraction
                                               └─ ratio columns
                                                    └─ CSV / figures / HTML report
```

---

## 2 · Function reference

### 2.1 `read_vsi(vsi_path)`

**Purpose:** Load the fluorescence series from a CellSens `.vsi` file.

**Parameters:**

| Name | Type | Description |
|---|---|---|
| `vsi_path` | `str` or `Path` | Absolute path to the `.vsi` file |

**Returns:** `(data, channel_names)` — `data` is a `(C, Y, X)` NumPy array, dtype `uint16`; `channel_names` is the list of channel name strings from BioFormats metadata.

**Internal mechanics:**

Uses `AICSImage` with the `BioformatsReader` backend (which wraps the Java BioFormats library via `bioformats_jar`). CellSens VSI files contain multiple series:

| Series | Contents | Shape | Action |
|---|---|---|---|
| S=0 | GFP + mCherry | 2 × 2048 × 2048 | **Loaded** |
| S=1 | Brightfield | 1 × 2048 × 2048 | Ignored |
| S=2 | Macro image | 3 × 512 × 512 | Ignored |

The call `img.get_image_data('CYX', S=0, T=0, Z=0)` selects Series 0 and returns a 3-D array with dimension order C (channel), Y, X. An assertion immediately verifies that exactly 2 channels are present — if the VSI has a non-standard structure the pipeline fails loudly rather than silently misassigning channels.

**Called in:** Section 2 (single-file load), Section 4c batch loop, Section 6 headless batch loop.

---

### 2.2 `split_channels(data, channel_names, gfp_idx=0, mcherry_idx=1)`

**Purpose:** Index into the CYX array to extract the two fluorescence channels.

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `data` | ndarray (C, Y, X) | — | Full multi-channel array from `read_vsi()` |
| `channel_names` | list of str | — | BioFormats channel names |
| `gfp_idx` | int | 0 | Index of the GFP plane within axis 0 |
| `mcherry_idx` | int | 1 | Index of the mCherry plane within axis 0 |

**Returns:** `(gfp_img, mcherry_img)` — two 2-D arrays `(Y, X)`, dtype `uint16`.

**Internal mechanics:** Simple NumPy indexing (`data[gfp_idx]`, `data[mcherry_idx]`). The function prints the full channel listing on every call so accidental swaps are immediately visible. Indices are user-configurable via `GFP_CHANNEL_IDX` / `MCHERRY_CHANNEL_IDX` in Section 0b.

**Called in:** Section 2, Section 4c loop, Section 6 loop.

---

### 2.3 `normalise(img)`

**Purpose:** Map raw fluorescence counts to the 0–255 float32 range using percentile stretching.

**Parameters:**

| Name | Type | Description |
|---|---|---|
| `img` | ndarray (Y, X) | Raw uint16 fluorescence image |

**Returns:** ndarray (Y, X), dtype `float32`, range [0, 255].

**Mathematics:**

Let $I$ denote the input image. Define:

$$\ell = P_1(I), \quad h = P_{99}(I)$$

where $P_k$ denotes the $k$-th percentile computed over all pixels. The output is:

$$\hat{I} = \text{clip}\!\left(\frac{I - \ell}{h - \ell + \varepsilon},\, 0,\, 1\right) \times 255$$

with $\varepsilon = 10^{-9}$ inserted to prevent division by zero when the image is uniform.

**Why percentile stretching rather than min–max?** Fluorescence images routinely contain hot pixels (saturated detectors) and dead pixels (zero counts). Using $P_1$ and $P_{99}$ instead of the true minimum and maximum makes the normalisation robust to these outliers. Clipping after scaling removes any pixel below $P_1$ or above $P_{99}$.

**Usage context:** Cellpose expects an input image whose intensities convey relative brightness; it does not expect raw photon counts. Normalise is also used for all matplotlib display panels in the diagnostic figures.

**Called in:** `segment_guvs()`, Section 2 channel sanity plot, Section 4c diagnostic figures, Section 5-iv evaluation.

---

### 2.4 `segment_guvs(gfp_img, diameter, flow_threshold, cellprob_threshold, model_type, use_gpu, min_size)`

**Purpose:** Run Cellpose instance segmentation on the normalised GFP channel and return an integer label mask.

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `gfp_img` | ndarray (Y, X) uint16 | — | Raw GFP image (normalised internally) |
| `diameter` | float or None | None | Expected object diameter in px. `None` = Cellpose auto-estimates |
| `flow_threshold` | float | 0.4 | Flow consistency post-filter threshold |
| `cellprob_threshold` | float | -2.0 | Foreground probability threshold |
| `model_type` | str | `'cyto3'` | Pretrained model identifier |
| `use_gpu` | bool | False | Enable CUDA acceleration |
| `min_size` | int | 15 | Minimum mask area in px²; smaller components are discarded |

**Returns:** `int32` ndarray (Y, X). Pixel value 0 = background; values 1…N are GUV labels.

**Cellpose algorithm — detailed:**

Cellpose uses a U-Net convolutional neural network. Its two-headed output jointly predicts:

1. **Cell-probability map** $p(x, y) \in [0, 1]$: sigmoid output estimating the probability that pixel $(x,y)$ belongs to any cell.

2. **Flow field** $\mathbf{F}(x, y) = (F_x, F_y)$: a 2-D vector at each pixel pointing in the direction of steepest ascent of the distance transform of the nearest cell mask. During training, ground-truth flows are derived analytically from labelled masks using the normalised gradient of a simulated heat diffusion from each cell interior.

**Foreground detection:**  
Pixels with $p(x, y) > \theta_{\text{prob}}$ (controlled by `cellprob_threshold`) are classified as foreground. Lower $\theta_{\text{prob}}$ → larger foreground region → detects dimmer objects.

- `cyto3` convention: $\theta_{\text{prob}} \approx -2$ (logit scale, mapped to probability $\approx 0.12$)
- `cpsam` convention: $\theta_{\text{prob}} = 0$ (logit scale, mapped to probability $= 0.5$)

**Mask assembly (gradient following):**  
For each foreground pixel, a particle is initialised at that position and advected along the predicted flow field:

$$\mathbf{x}(t+1) = \mathbf{x}(t) + \mathbf{F}(\mathbf{x}(t))$$

Particles whose trajectories converge to the same spatial attractor are assigned the same label. This is iterated for a fixed number of steps (typically 200). The result is a label image where connected groups of pixels sharing an attractor form one instance mask.

**Flow consistency filter (`flow_threshold`):**  
After mask assembly, Cellpose recomputes the flow field implied by the assembled mask (by diffusion from the mask pixels) and compares it to the network's predicted flow:

$$e_{\mathrm{flow}} = \left\|\mathbf{F}_{\mathrm{pred}} - \mathbf{F}_{\mathrm{mask}}\right\|^2$$

Masks whose mean $e_{\mathrm{flow}}$ exceeds `flow_threshold` are discarded as inconsistent. Lower `flow_threshold` → more permissive → retains masks with noisier flow predictions.

**`min_size` filter:**  
After all above steps, any connected component with area smaller than `min_size` px² is removed. In this pipeline `min_size = CP_MIN_AREA_PX` which is derived from the physical minimum diameter:

$$d_{\min}^{\text{px}} = d_{\min}^{\mu\text{m}} \times R_{\text{px}/\mu\text{m}} = 8.0 \times 15.3846 \approx 123\,\text{px}$$

$$A_{\min} = \pi \left(\frac{d_{\min}^{\text{px}}}{2}\right)^2 = \pi \times 61.5^2 \approx 11{,}882\,\text{px}^2$$

**`diameter` parameter (cyto3 only):**  
Cellpose internally rescales the input image so that objects appear at a canonical size (~30 px diameter in the network's receptive field). When `diameter=None`, Cellpose auto-estimates the scale factor from the image. When `diameter` is supplied it directly sets the rescale factor. cpsam is size-invariant and ignores this parameter.

**Mean diameter log (when `diameter=None`):**  
After segmentation the pipeline computes the equivalent circular diameter of each detected mask:

$$\hat{d}_i = 2\sqrt{A_i / \pi}$$

and prints the mean $\bar{d} = \frac{1}{N}\sum_{i=1}^{N} \hat{d}_i$ so the user can set `CP_DIAMETER` for future runs.

**Section 5-iii redefinition:**  
When `USE_FINETUNED_MODEL = True`, Section 5-iii redefines `segment_guvs()` in-place to load the fine-tuned weights. All downstream code (Section 4c, Section 6) automatically uses the new implementation without changes.

**Called in:** Section 4c loop (interactive), Section 6 loop (headless). Also called inside the archived Section 3a raw cell via `segment_guvs_sam()`.

---

### 2.5 `make_membrane_lumen_masks(masks, membrane_thickness_px=5)`

**Purpose:** Decompose each GUV label into a membrane ring and a lumen interior by repeated binary erosion.

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `masks` | ndarray (Y, X) int32 | — | Cellpose label mask |
| `membrane_thickness_px` | int | 5 | Number of erosion iterations |

**Returns:** Two dicts `{label (int) → bool ndarray (Y, X)}`: `membrane_masks` and `lumen_masks`.

**Mathematics — binary morphological erosion:**

For label $\ell$, let $M$ denote the integer label mask array (0 = background). Define the binary mask:

$$B_\ell(x,y) = \begin{cases} 1 & \text{if } M(x,y) = \ell \\ 0 & \text{otherwise} \end{cases}$$

The structuring element $S$ is the 4-connected cross (generated by `ndi.generate_binary_structure(2, 1)`):

$$S = \begin{pmatrix} 0 & 1 & 0 \\ 1 & 1 & 1 \\ 0 & 1 & 0 \end{pmatrix}$$

The single-step erosion of $B$ by $S$ is:

$$(B \ominus S)(x,y) = 1 \iff S_{(x,y)} \subseteq B$$

i.e., a pixel survives erosion if and only if every pixel covered by $S$ centred at that pixel is also in $B$. Equivalently, a pixel is removed if any of its 4-neighbours (up, down, left, right) lies outside $B$.

After $t$ iterations (where $t$ = `membrane_thickness_px`):

$$B_\ell^{(t)} = B_\ell \ominus S \ominus S \ominus \cdots \quad (t \text{ times})$$

**Membrane ring:**

$$\text{membrane}_\ell = B_\ell \setminus B_\ell^{(t)} = B_\ell \,\text{ AND NOT }\, B_\ell^{(t)}$$

This is the annular region stripped away by erosion — the outer shell of each GUV.

**Lumen:**

$$\text{lumen}_\ell = B_\ell^{(t)}$$

The surviving interior after $t$ erosion steps.

**Geometric interpretation:**  
Erosion by $t$ steps of a 4-connected cross is equivalent to removing all pixels within a Chebyshev distance of $t$ from the mask boundary (in the $L^1$ metric under 4-connectivity). For a circular GUV of radius $r$ px, the lumen radius is approximately $r - t$ px. The membrane ring has approximate width $t$ px.

**Edge case — GUV too small for full erosion:**  
If $B_\ell^{(t)}$ is empty (all pixels removed, i.e., $r < t$), then:
- `membrane_masks[label] = B_ℓ` (full mask)
- `lumen_masks[label] = np.zeros_like(B_ℓ)` (empty)

This prevents a degenerate state where no region is quantifiable. The resulting `GFP_mem_lumen_ratio` will be `NaN` because `l_px.size == 0` triggers the guard in `quantify()`.

**Called in:** `run_quantification()`, Section 4c ring-preview figure, Section 4b archived raw cell.

---

### 2.6 `quantify(image, masks, membrane_masks, lumen_masks, channel_name)`

**Purpose:** Measure intensity statistics in the membrane and lumen of every GUV for one channel.

**Parameters:**

| Name | Type | Description |
|---|---|---|
| `image` | ndarray (Y, X) uint16 | Raw fluorescence image for one channel |
| `masks` | ndarray (Y, X) int32 | Full label mask |
| `membrane_masks` | dict {int → bool ndarray} | From `make_membrane_lumen_masks()` |
| `lumen_masks` | dict {int → bool ndarray} | From `make_membrane_lumen_masks()` |
| `channel_name` | str | Column name prefix (`'GFP'` or `'mCherry'`) |

**Returns:** `pd.DataFrame` with one row per GUV label.

**Mathematics:**

For label $\ell$, let:
- $\Omega_\ell^m$ = set of pixels in the membrane mask
- $\Omega_\ell^l$ = set of pixels in the lumen mask
- $\Omega_\ell$ = full GUV mask ($\Omega_\ell^m \cup \Omega_\ell^l$)
- $I(x,y)$ = raw pixel value (uint16 photon count)

**Centroid:**

$$\bar{x}_\ell = \frac{1}{|\Omega_\ell|} \sum_{(x,y)\in\Omega_\ell} x, \qquad \bar{y}_\ell = \frac{1}{|\Omega_\ell|} \sum_{(x,y)\in\Omega_\ell} y$$

**Area and equivalent diameter:**

$$A_\ell = |\Omega_\ell|, \qquad d_\ell = 2\sqrt{A_\ell / \pi}$$

The diameter formula inverts the circle-area relationship $A = \pi r^2$, giving the diameter of a circle with the same pixel count. This is the **equivalent circular diameter** — a size metric that does not assume circularity but is expressed in familiar units.

**Mean intensity:**

$$\mu_m = \frac{1}{|\Omega_\ell^m|} \sum_{(x,y)\in\Omega_\ell^m} I(x,y), \qquad \mu_l = \frac{1}{|\Omega_\ell^l|} \sum_{(x,y)\in\Omega_\ell^l} I(x,y)$$

Mean intensity is **extensive with pixel area**: it represents the mean photon count per pixel in the region, directly comparable across GUVs of different sizes.

**Integrated intensity:**

$$\Sigma_m = \sum_{(x,y)\in\Omega_\ell^m} I(x,y) = |\Omega_\ell^m| \cdot \mu_m$$

Integrated (total) intensity is proportional to the total number of fluorophores in the region. It is **not** corrected for region size, so it grows with GUV area and is best used for within-GUV comparisons or when normalised explicitly.

**Membrane-to-lumen ratio:**

$$R_\ell = \frac{\mu_m}{\mu_l + \varepsilon}$$

with $\varepsilon = 10^{-9}$. The ratio quantifies how much more concentrated the fluorophore is at the membrane relative to the lumen.

- $R > 1$: membrane-enriched (protein or lipid localises to the membrane)
- $R = 1$: uniform distribution
- $R < 1$: lumen-enriched

The ratio is set to `NaN` when:
- `l_px.size == 0` (lumen region is empty — GUV too small to erode)
- `denom <= 0` (lumen mean is zero or negative — pathological case)

The $\varepsilon$ term only prevents numerical infinity in the edge case where the lumen mean is non-zero but extremely small; it does not replace the `NaN` guards above.

**Called in:** `run_quantification()` (twice — once per channel).

---

### 2.7 `run_quantification(gfp_img, mcherry_img, masks, membrane_thickness_px=5)`

**Purpose:** Full quantification: from a label mask produce a merged per-GUV DataFrame for both channels.

**Parameters:**

| Name | Type | Description |
|---|---|---|
| `gfp_img` | ndarray (Y, X) uint16 | GFP channel |
| `mcherry_img` | ndarray (Y, X) uint16 | mCherry channel |
| `masks` | ndarray (Y, X) int32 | Cellpose/corrected label mask |
| `membrane_thickness_px` | int | Passed to `make_membrane_lumen_masks()` |

**Returns:** `(df, membrane_masks, lumen_masks)` where `df` is the merged DataFrame.

**Internal logic:**

1. Calls `make_membrane_lumen_masks()` once (shared rings for both channels).
2. Calls `quantify()` for GFP → `df_gfp`.
3. Calls `quantify()` for mCherry → `df_mcherry`.
4. Merges on the `label` column. Geometric columns (`label`, `centroid_x_px`, `centroid_y_px`, `area_px2`, `diameter_px`) come from `df_gfp`; channel-specific intensity columns come from both frames.

The merge guarantees that both channels report identical geometry columns (they are computed from the same `masks` array) and that each row represents one GUV with all measurements present.

**Called in:** Section 4c batch loop, Section 6 headless loop.

---

### 2.8 `_wait_for_close(viewer)` *(Section 4c)*

**Purpose:** Block Python execution until the napari window is closed, by pumping the Qt event loop.

**Parameters:**

| Name | Type | Description |
|---|---|---|
| `viewer` | `napari.Viewer` | Active napari viewer instance |

**Returns:** None.

**Internal mechanics:**

```python
_app  = QApplication.instance()
_qwin = viewer.window._qt_window
while _qwin.isVisible():
    _app.processEvents()
    time.sleep(0.05)
```

Inside a Jupyter kernel, the Qt event loop is not running autonomously. Without event-loop processing, the napari window would freeze (non-responsive to mouse/keyboard). The while loop polls `isVisible()` every 50 ms, calling `processEvents()` each iteration to allow Qt to handle GUI events (repaints, clicks, keyboard input). Once the user closes the window, `isVisible()` returns `False` and control returns to the Python caller.

**Called in:** Section 4c batch loop (once per file, after opening the napari viewer).

---

### 2.9 `_bulk_mode(img, mask_ok, bins=512)` *(inline function, Sections 4c and 6)*

**Purpose:** Estimate the background fluorescence level of the external aqueous phase by taking the mode of the out-of-mask pixel intensity histogram.

**Parameters:**

| Name | Type | Default | Description |
|---|---|---|---|
| `img` | ndarray (Y, X) uint16 | — | Raw fluorescence image (GFP or mCherry) |
| `mask_ok` | ndarray (Y, X) int32 | — | Combined GUV label mask; pixels with value 0 are outside all GUVs |
| `bins` | int | 512 | Number of histogram bins |

**Returns:** `float` — estimated background intensity (mode of the out-of-mask distribution).

**Mathematics:**

**Step 1 — Extract background pixels:**

$$\mathcal{B} = \bigl\{\, I(x,y) : M(x,y) = 0 \,\bigr\}$$

where $M$ is the label mask (0 outside all GUVs). These are all pixels belonging to the external aqueous phase (plus any field-of-view area not covered by a GUV mask).

**Step 2 — Histogram:**

The range $[I_{\min}(\mathcal{B}),\, I_{\max}(\mathcal{B})]$ is divided into 512 equal-width bins. Let $c_k$ be the count in bin $k$, and $[e_k, e_{k+1})$ be the bin edge interval.

$$k^{\ast} = \underset{k}{\arg\max}\; c_k$$

**Step 3 — Bin-centred mode estimate:**

$$\hat{\mu}_{\text{bg}} = \frac{e_{k^{\ast}} + e_{k^{\ast}+1}}{2}$$

The midpoint of the peak bin is taken as the mode estimate. With 512 bins over a uint16 range ($[0, 65535]$), each bin spans approximately 128 intensity units, giving a background estimate accurate to ±64 counts.

**Statistical rationale — why the mode?**

The external aqueous phase typically occupies the majority of the image area (GUVs are sparse). The out-of-mask pixel distribution is therefore strongly peaked at the autofluorescence + detector offset baseline and has a long right tail from any debris, vesicle boundaries, or bright background features. The **mode** (most frequent value) locates this peak robustly and is not pulled toward the right tail the way the **mean** would be. The **median** is also robust but requires sorted data; the histogram approximation is faster on large pixel arrays.

**Called in:** Section 4c loop (after each napari correction session, applied to both GFP and mCherry), Section 6 headless loop (same position).

---

### 2.10 `segment_guvs_sam(gfp_img, flow_threshold, cellprob_threshold, use_gpu, min_size)` *(archived Section 3a)*

An archived variant of `segment_guvs()` specifically for the `cpsam` model. It omits the `diameter` and `channels` parameters because Cellpose-SAM is size-invariant and channel-order invariant. This function is kept in a Raw cell for reference; the main pipeline uses `segment_guvs()` with `model_type='cpsam'` in Section 4c/6 instead.

---

## 3 · Background subtraction and dual-channel ratio

After calling `run_quantification()`, Sections 4c and 6 compute three additional derived columns per GUV using `_bulk_mode()`:

**Background-subtracted lumen intensities:**

$$I_{\text{GFP,lumen}}^{\text{bgsub}} = \max\!\left(0,\; \mu_l^{\text{GFP}} - \hat{\mu}_{\text{bg}}^{\text{GFP}}\right)$$

$$I_{\text{mCh,lumen}}^{\text{bgsub}} = \max\!\left(0,\; \mu_l^{\text{mCherry}} - \hat{\mu}_{\text{bg}}^{\text{mCherry}}\right)$$

The `.clip(0)` (equivalent to $\max(0, \cdot)$) prevents negative values. A negative value arises when the lumen mean is below the background level, which is consistent with zero net lumen fluorescence above baseline.

**Dual-channel lumen ratio:**

$$R_{\text{lumen}} = \frac{I_{\text{GFP,lumen}}^{\text{bgsub}}}{I_{\text{mCh,lumen}}^{\text{bgsub}} + \varepsilon}$$

with $\varepsilon = 10^{-9}$.

**Why background subtraction is necessary for this ratio:**

Without background subtraction, the ratio $\mu_l^{\text{GFP}} / \mu_l^{\text{mCherry}}$ would be confounded by differential autofluorescence and detector offsets between the two channels. The GFP and mCherry detectors have different quantum efficiencies, gain settings, and dark-noise levels. Subtracting the channel-specific background mode removes the additive baseline from each channel independently before dividing, making the ratio a measure of the relative lumen concentrations of the two labelled species above baseline.

The ratio is stoichiometrically interpretable as:

$$R_{\text{lumen}} \approx \frac{[\text{GFP-labelled species}]_{\text{lumen}}}{[\text{mCherry-labelled species}]_{\text{lumen}}}$$

provided that both fluorophores have comparable extinction coefficients and quantum yields under the imaging conditions (or that a calibration factor is applied externally).

**Why the lumen rather than the membrane for this ratio:**

The lumen encapsulates the soluble interior content of the GUV. The lumen ratio directly probes the relative concentrations of soluble cargo. The membrane ratio (`GFP_mem_lumen_ratio`, `mCherry_mem_lumen_ratio`) instead measures how much of each channel's signal is enriched at the membrane relative to the interior, which is a probe of membrane localisation.

---

## 4 · Statistical and analytical methods

### 4.1 Membrane / lumen ratio as a localisation index

The ratio $R = \mu_{\text{mem}} / (\mu_{\text{lumen}} + \varepsilon)$ is a **dimensionless enrichment index**:

- $R \gg 1$: fluorophore is membrane-localised (e.g., a lipid-anchored protein)
- $R \approx 1$: fluorophore is uniformly distributed
- $R \ll 1$: fluorophore is excluded from the membrane (e.g., a large cytosolic protein)

The reference line at $R = 1$ in the scatter plot (Section 4c summary figure) marks this boundary. The quadrant annotations count GUVs in each quadrant of the $(R_{\text{GFP}}, R_{\text{mCherry}})$ plane:

| Quadrant | Condition | Biological interpretation |
|---|---|---|
| Top-right | $R_{\text{GFP}} > 1$ AND $R_{\text{mCh}} > 1$ | Both species membrane-enriched |
| Bottom-right | $R_{\text{GFP}} > 1$ AND $R_{\text{mCh}} < 1$ | GFP membrane-enriched, mCherry lumen-enriched |
| Top-left | $R_{\text{GFP}} < 1$ AND $R_{\text{mCh}} > 1$ | Opposite pattern |
| Bottom-left | $R_{\text{GFP}} < 1$ AND $R_{\text{mCh}} < 1$ | Both lumen-enriched |

### 4.2 Violin plot distributions (Section 4c summary figure)

The seaborn violin plot uses **Kernel Density Estimation (KDE)** to estimate the probability density of ratio values:

$$\hat{f}(x) = \frac{1}{Nh} \sum_{i=1}^{N} K\!\left(\frac{x - x_i}{h}\right)$$

where $K$ is a Gaussian kernel and $h$ is the bandwidth (Scott's rule by default). The `cut=0` argument truncates the KDE at the observed data range (prevents density from extending beyond real measurements). The `inner='box'` overlay draws a standard box plot inside the violin:

- Central line: median ($Q_2$)
- Box: interquartile range $[Q_1, Q_3]$
- Whiskers: $Q_1 - 1.5 \cdot \text{IQR}$ to $Q_3 + 1.5 \cdot \text{IQR}$

### 4.3 Diameter vs ratio scatter (Section 4c summary figure)

This plot visualises whether the membrane/lumen ratio depends on GUV size. A significant correlation would indicate a systematic artefact, most likely from the fixed `MEMBRANE_THICKNESS_PX` parameter:

- For small GUVs: a fixed membrane thickness represents a proportionally larger fraction of the total GUV volume, causing the lumen to shrink disproportionately and the ratio to be upward-biased.
- For large GUVs: the membrane ring is thin relative to the total area, and the ratio is less sensitive to `MEMBRANE_THICKNESS_PX`.

If a strong size-dependence is observed, consider scaling `MEMBRANE_THICKNESS_PX` proportionally to GUV diameter, or restricting analysis to a narrower size bin.

### 4.4 Average Precision @ IoU = 0.5 (Section 5-iv evaluation)

This metric evaluates the quality of the fine-tuned model against a ground-truth mask on a held-out image.

**Intersection over Union (IoU):**

For a predicted mask $\hat{M}_i$ and ground-truth mask $G_j$:

$$\text{IoU}(\hat{M}_i, G_j) = \frac{|\hat{M}_i \cap G_j|}{|\hat{M}_i \cup G_j|}$$

IoU = 1 for a perfect overlap; IoU = 0 for no overlap. The threshold of 0.5 is the standard COCO detection criterion ("50% overlap").

**Matching:**

At threshold $\tau = 0.5$, each predicted mask is matched to at most one ground-truth mask (the one with the highest IoU, provided IoU $\geq \tau$). The matching is done greedily in descending order of prediction confidence.

- **True Positive (TP):** predicted mask matched to a ground truth (IoU $\geq 0.5$)
- **False Positive (FP):** predicted mask with no matching ground truth
- **False Negative (FN):** ground-truth mask with no matching prediction

**Precision and Recall:**

$$\text{Precision} = \frac{\text{TP}}{\text{TP} + \text{FP}}, \qquad \text{Recall} = \frac{\text{TP}}{\text{TP} + \text{FN}}$$

**Average Precision:**

$$\text{AP} = \sum_k \bigl(\text{Recall}(k) - \text{Recall}(k-1)\bigr) \cdot \text{Precision}(k)$$

This is the area under the precision–recall curve, computed at the single threshold $\tau = 0.5$ here (not averaged across thresholds as in the full COCO metric). AP = 1.0 means every GUV was correctly detected with IoU $\geq 0.5$; AP = 0.0 means no correct detections.

### 4.5 Cellpose fine-tuning (Section 5-ii)

Fine-tuning adapts the pretrained Cellpose weights to GUV-specific morphology using the manually corrected masks as ground truth.

**Optimiser:** Stochastic Gradient Descent (SGD).

**Loss function:**  
Cellpose training minimises a composite loss combining cell-probability cross-entropy and flow-field mean squared error:

$$\mathcal{L}_{\text{total}} = \mathcal{L}_{\text{prob}} + \lambda_{\text{flow}} \mathcal{L}_{\text{flow}}$$

$$\mathcal{L}_{\text{prob}} = -\frac{1}{N}\sum_{(x,y)} \left[ y\log p(x,y) + (1-y)\log(1-p(x,y)) \right]$$

$$\mathcal{L}_{\text{flow}} = \frac{1}{N}\sum_{(x,y)\in\text{foreground}} \left\| \hat{\mathbf{F}}(x,y) - \mathbf{F}^*(x,y) \right\|^2$$

where $y$ is the ground-truth foreground label (1 inside cells, 0 outside), $p$ is the predicted cell probability, $\hat{\mathbf{F}}$ is the predicted flow field, and $\mathbf{F}^*$ is the ground-truth flow derived from the annotated masks by heat-diffusion simulation.

**Weight update rule (SGD with L2 regularisation):**

$$w_j^{(t+1)} = w_j^{(t)} - \eta \left[ \frac{\partial \mathcal{L}_{\text{total}}}{\partial w_j} + \lambda w_j^{(t)} \right]$$

Equivalently (factoring the weight term):

$$w_j^{(t+1)} = w_j^{(t)} (1 - \eta \lambda) - \eta \frac{\partial \mathcal{L}_{\text{total}}}{\partial w_j}$$

**Parameters:**

| Parameter | Value | Role |
|---|---|---|
| $\eta$ (learning rate) | 0.1 | Step size per gradient update |
| $\lambda$ (weight decay) | $10^{-5}$ | L2 penalty coefficient |
| Epochs | 100 | Full passes over the training set |

**L2 regularisation (weight decay):**  
The term $\lambda w_j^{(t)}$ in the gradient shrinks each weight toward zero at every update. This penalises large-magnitude weights, discouraging the model from memorising the training images and encouraging smoother, more generalisable decision boundaries. At $\lambda = 10^{-5}$ the penalty is mild — it primarily prevents unbounded weight growth rather than aggressively constraining the model.

---

## 5 · Output figure descriptions

### Section 4c per-file figures

**`{stem}_segmentation_preview.png` — 4 panels:**
1. GFP normalised (Greens colormap)
2. mCherry normalised (Reds colormap)
3. Label mask coloured by `nipy_spectral` (one colour per GUV), with integer label IDs overlaid at each centroid
4. RGB composite: red channel = mCherry normalised, green channel = GFP normalised, blue = 0; mask overlaid at 30% opacity

**`{stem}_ring_preview.png` — 4 panels:**
1. Label mask (nipy_spectral)
2. Ring annotation on black background: orange = membrane pixels, blue = lumen pixels
3. GFP (Greens) + ring annotation at 55% opacity
4. mCherry (Reds) + ring annotation at 55% opacity

The RGBA overlay encodes per-pixel region class rather than per-GUV colour. This allows visual verification that the orange ring aligns with the visible fluorescence ring in the GFP channel.

**`{stem}_ratio_scatter.png`:**  
Scatter of $R_{\text{GFP}}$ (x-axis) vs $R_{\text{mCherry}}$ (y-axis) per GUV within the file. Useful for identifying the distribution of localisation phenotypes within a single image.

**`{stem}_lumen_ratio.png`:**  
Scatter of `mCherry_lumen_bgsub` (x-axis) vs `GFP_lumen_bgsub` (y-axis) per GUV, with point colour encoding `GFP_mCherry_lumen_ratio` on a RdBu_r diverging colormap (centre = 1, range 0–2). A dashed diagonal marks the 1:1 line. Points above the diagonal have more GFP than mCherry in the lumen; points below have more mCherry.

### Section 4c combined figure

**`batch_interactive_lumen_ratio.png` — 3 panels:**
1. $R_{\text{GFP}}$ vs $R_{\text{mCherry}}$ scatter coloured by file, with quadrant count annotations and reference lines at ratio = 1
2. Violin + box plot of ratio distributions, one violin per channel (GFP and mCherry), pooled across all files
3. GUV diameter (px) vs $R_{\text{mCherry}}$ scatter coloured by file, with reference line at ratio = 1

---

## 6 · Data flow and variable lifetimes

| Variable | Created in | Consumed in | Type |
|---|---|---|---|
| `vsi_files` | Section 0b | Sections 2, 4c, 6 | `list[Path]` |
| `gfp_img`, `mcherry_img` | Section 2 | Sections 4 (archived), 4b (archived) | `ndarray (Y, X) uint16` |
| `masks` | Section 4 (archived) | Section 4 (archived) | `ndarray (Y, X) int32` |
| `CHOSEN_MODEL` | Section 3a | Sections 4c, 6 | `str` |
| `batch_df` | Section 4c | Section 7 | `pd.DataFrame` |
| `_training_pairs` | Section 4c loop | Section 5-i | `list[(stem, gfp_array, masks_array)]` |
| `_report_figs` | Section 4c loop | Section 4c HTML generation | `list[(title_str, Path)]` |
| `all_dfs` | Section 6 | Section 7 | `list[pd.DataFrame]` |
| `combined` | Section 7 | User (final CSV) | `pd.DataFrame` |
| `FINETUNED_MODEL_PATH` | Section 5-ii | Sections 5-iii, 5-iv | `Path` |

---

## 7 · Key design decisions and their rationale

**1. Segmentation channel = GFP (configurable).**  
GFP is typically expressed as a membrane-targeted construct in GUV experiments and produces clear ring signal. The GFP ring has higher signal-to-noise than mCherry in many setups because GFP has a higher quantum yield (~0.6 vs ~0.22 for mCherry). The segmentation channel is configurable if mCherry provides a cleaner boundary.

**2. Erosion with a 4-connected structuring element (cross), not 8-connected (square).**  
The 4-connected element removes only pixels where a direct (axis-aligned) neighbour is outside the mask. This is more conservative than the 8-connected (square) element and produces a membrane ring that more closely tracks the actual mask boundary contour at diagonal edges. It avoids removing corner pixels that would otherwise be retained under an 8-connected criterion.

**3. Background mode from out-of-mask pixels, not a blank region.**  
Computing background from a fixed blank region of interest assumes that the external phase has uniform intensity everywhere in the image. In practice, there can be spatial gradients from optical non-uniformity. Using all out-of-mask pixels averages over the entire field, giving a more representative estimate of the mean background. For structured backgrounds, a spatially local estimate (e.g., a local histogram per GUV centroid) would be more accurate but is not implemented here.

**4. $\varepsilon = 10^{-9}$ as denominator stabiliser, not a biological floor.**  
The small constant prevents division by zero in the ratio computation but does not represent any physical minimum intensity. Ratios computed when the denominator is close to $\varepsilon$ are biologically meaningless and should be treated as effectively infinite. The explicit `NaN` guard (`l_px.size == 0`) catches the dominant case (empty lumen) before the division.

**5. `_bulk_mode` is redefined identically in Sections 4c and 6.**  
The function is an inline closure (not a module-level function) because it is only three lines long and is always called immediately after use. This avoids polluting the module namespace and makes the context of each call explicit.

**6. Section 5-iii redefines `segment_guvs()` in-place.**  
Rather than adding a conditional branch inside `segment_guvs()`, the function is entirely replaced. This keeps the call sites in Sections 4c and 6 unchanged and avoids an extra global flag lookup on every call.

---

## 8 · Limitations and known approximations

1. **Membrane ring has fixed physical thickness.** `MEMBRANE_THICKNESS_PX` is a constant regardless of GUV size. For very large GUVs, the ring represents a smaller fraction of GUV radius; for very small GUVs, it may consume the entire interior. The diameter vs ratio plot in the summary figure is designed to flag this artefact.

2. **Background estimated once per image.** A single global background value is subtracted from all GUVs in the image. If GUVs are densely packed or if there is a spatial gradient in background, this approximation introduces per-GUV errors.

3. **Mode estimated from binned histogram.** The bin-centred mode has a precision of ±(intensity range / 2 × 512) counts. For uint16 data (~65535 range) this is ±64 counts. For very low-signal images or very clean images with narrow intensity histograms, increasing `bins` (e.g., to 1024) improves precision.

4. **Equivalent circular diameter assumes circularity.** `diameter_px = 2√(area/π)` gives the diameter of a circle with the same area as the mask. Elongated or irregular GUVs will have a `diameter_px` that underestimates their longest axis.

5. **Integrated intensity not corrected for unequal mask areas.** `GFP_membrane_integrated` is the sum of pixel values in the membrane ring. The membrane ring area varies between GUVs (larger GUVs have proportionally larger membrane rings), so integrated intensity is not directly comparable across GUVs of different sizes. Use `_mean` columns for size-independent comparisons.

6. **No photobleaching correction.** If multiple files are imaged sequentially with long acquisition times, later files may show reduced fluorescence due to photobleaching. No time-ordering correction is applied.

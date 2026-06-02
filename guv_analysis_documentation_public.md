# GUV Fluorescence Analysis Pipeline — Documentation

**Script:** `guv_analysis.ipynb`  
**Language:** Python 3.10 (Jupyter Notebook)

---

## Overview

This pipeline segments and quantifies fluorescence from two-channel (GFP + mCherry) confocal images of Giant Unilamellar Vesicles (GUVs) stored as CellSens `.vsi` files. It uses [Cellpose](https://www.cellpose.org/) to segment GUV boundaries on the GFP channel, allows manual mask correction in [napari](https://napari.org/), then measures fluorescence intensity separately in the **membrane ring** and **lumen** of each GUV.

### Two workflows

| Workflow | Sections | Use case |
|---|---|---|
| **Interactive batch** | 0b → 1 → 2 → 3a → 4c | Manually correct every file in napari; recommended for new datasets |
| **Headless batch** | 0b → 1 → 2 → 3a → 4c (subset) → 5 → 6 → 7 | Tune on a few files interactively, then auto-process the rest |

---

## Installation

```bash
conda create -n guv_env python=3.10
conda activate guv_env
pip install aicsimageio[bioformats] bioformats_jar cellpose "napari[all]" pyqt5 \
            numpy scipy pandas matplotlib seaborn tifffile notebook ipykernel
python -m ipykernel install --user --name guv_env --display-name "GUV pipeline"
```

**Requirements:**
- Java JDK 8+ on PATH (for BioFormats VSI reading)
- macOS: Xcode Command Line Tools (`xcode-select --install`)
- GPU (optional): CUDA-compatible GPU for faster Cellpose inference

**Download Cellpose-SAM model on first use (~400 MB):**
```python
from cellpose import models
models.CellposeModel(pretrained_model='cpsam')
```

---

## Configuration (Section 0b)

All user-editable parameters live in a single cell. Edit these before running the pipeline.

### Paths

| Variable | Description |
|---|---|
| `INPUT_DIR` | Folder containing `.vsi` files |
| `OUTPUT_DIR` | Destination for all output files; created automatically if absent |

### Channel indices

| Variable | Default | Description |
|---|---|---|
| `GFP_CHANNEL_IDX` | `0` | Index of the GFP channel within Series 0 of the VSI file |
| `MCHERRY_CHANNEL_IDX` | `1` | Index of the mCherry channel |
| `SEGMENTATION_CHANNEL_IDX` | `GFP_CHANNEL_IDX` | Which channel to segment on; change to `MCHERRY_CHANNEL_IDX` if mCherry has cleaner ring signal or contain cell marker (i.e. dextran-TexasRed)|

Run the **channel sanity check** cell (Section 2) once per new dataset to visually confirm the assignments are correct. Flip the two values if the channels appear swapped.

### Cellpose segmentation parameters

| Variable | Default | Description |
|---|---|---|
| `SPATIAL_RESOLUTION_PX_UM` | `15.3846` | Image pixel size in px/µm; used to compute the minimum detection size |
| `MINIMUM_PHYSICAL_DIAMETER_UM` | `8.0` | Smallest GUV diameter to detect, in µm |
| `CP_MIN_DIAMETER_PX` | *(computed)* | `MINIMUM_PHYSICAL_DIAMETER_UM × SPATIAL_RESOLUTION_PX_UM` |
| `CP_MIN_AREA_PX` | *(computed)* | Minimum GUV area in px²; passed as `min_size` to Cellpose |
| `CP_DIAMETER` | `None` | Expected GUV diameter in pixels. `None` = Cellpose auto-estimates; set to the printed mean diameter to speed up subsequent runs |
| `CP_FLOW_THRESHOLD` | `0` | Cellpose flow error threshold. Lower → keep more detections |
| `CP_CELLPROB_THRESHOLD` | `0` | Cellpose cell probability threshold. Lower → detect dimmer/smaller objects |
| `CP_MODEL` | `'cyto3'` | Fallback model string used when `CHOSEN_MODEL != 'cpsam'` |
| `CP_GPU` | `True` | Enable GPU inference. Set `False` if no CUDA GPU is available |

### Quantification

| Variable | Default | Description |
|---|---|---|
| `MEMBRANE_THICKNESS_PX` | `7` | Thickness of the membrane ring eroded inward from each GUV boundary. Verify visually in Section 4b before running the batch |

---

## Workflow

### Section 1 · Import libraries and define helper functions

Loads all Python imports and defines the six core helper functions (`read_vsi`, `split_channels`, `normalise`, `segment_guvs`, `make_membrane_lumen_masks`, `run_quantification`). Must be run before any other section.

### Section 2 · Load a single file

**User variable:** `FILE_INDEX` — integer index into the sorted list of `.vsi` files found in `INPUT_DIR` (0-based).

Loads the file with `read_vsi()` and splits channels with `split_channels()`, producing `gfp_img` and `mcherry_img` (2-D NumPy arrays, dtype uint16). Two diagnostic plots are produced:

- **Channel sanity check** — greyscale panels for each channel so you can confirm GFP vs mCherry assignment
- **Quick colour preview** — GFP in green, mCherry in red

### Section 3a · Model comparison: cyto3 vs Cellpose-SAM *(archived)*

This cell is kept as a **Raw cell** (not executed by default). Convert it back to Code if you want to run a side-by-side napari comparison of the two models on the current file.

After comparison, set the model choice in the next cell:

```python
CHOSEN_MODEL = 'cpsam'   # or 'cyto3'
```

This variable is read by the batch loop in Section 4c and the headless loop in Section 6.

---

### Section 4 · Inspect & manually correct masks in napari *(archived)*

The individual single-file napari correction cells (4, 4a, 4b) are kept as **Raw cells** for reference. The recommended workflow is **Section 4c**, which wraps the same steps in a batch loop.

---

### Section 4c · Batch interactive processing

**User variable:** `FILE_INDICES` — list of integer indices into `vsi_files`.

```python
FILE_INDICES = [0, 1, 2]   # process these three files
```

Negative indices work (`-1` = last file, etc.). The loop prints which files will be processed before starting; adjust `FILE_INDICES` and re-run the setup cell if needed.

For each file the loop:
1. Loads and splits channels
2. Segments with `CHOSEN_MODEL` and the config parameters
3. Opens a napari window — edit the **Masks** label layer, then **close the window** to continue to the next file
4. Generates four figures (saved as PNG to `OUTPUT_DIR`):
   - `{stem}_segmentation_preview.png` — 4-panel: GFP | mCherry | masks with IDs | GFP+mCherry+masks composite overlay
   - `{stem}_ring_preview.png` — 4-panel: masks | membrane/lumen ring annotation | GFP overlay | mCherry overlay
   - `{stem}_ratio_scatter.png` — GFP vs mCherry membrane/lumen ratio per GUV
   - `{stem}_lumen_ratio.png` — background-subtracted lumen GFP vs mCherry scatter (colour = ratio)
5. Saves per-file measurements to `{stem}_measurements.csv`

After the loop, two additional outputs are produced:

- `batch_interactive_measurements.csv` — all files concatenated
- `batch_interactive_lumen_ratio.png` — three-panel combined figure: membrane/lumen ratio scatter (both channels) | violin distributions | mCherry ratio vs GUV diameter
- `batch_interactive_report.html` — self-contained HTML report with all figures embedded as base64 images, grouped by file

**Napari editing cheatsheet:**

| Action | How |
|---|---|
| Select a GUV | Alt-click on it, or type its label ID in the Labels panel |
| Paint a new GUV | Pick an unused label number, draw with brush |
| Erase a mask | Switch to label `0`, paint over it |
| Fill a closed outline | Bucket (fill) tool |
| Undo | Ctrl+Z |

The label number shown in the spinner is the only thing that determines which GUV you're drawing — napari will not prevent painting one GUV over another, so always pick the correct label first.

---

### Section 5 · Fine-tune Cellpose *(optional)*

Use napari-corrected masks as training data to adapt Cellpose to your specific GUVs. Run after accumulating masks from multiple correction sessions. Aim for ≥ 20 image/mask pairs before fine-tuning.

#### 5-i · Save training pairs

Saves all napari-corrected image/mask pairs accumulated during the Section 4c loop as TIFF files to `OUTPUT_DIR/training_data/`. Run once after Section 4c completes.

#### 5-ii · Fine-tune

| Variable | Default | Description |
|---|---|---|
| `START_MODEL` | `'cyto3'` | Base model to fine-tune from (`'cyto3'` or `'cpsam'`) |
| `N_EPOCHS` | `100` | Training epochs |
| `LEARNING_RATE` | `0.1` | SGD learning rate |
| `WEIGHT_DECAY` | `1e-5` | L2 regularisation |
| `FINETUNED_MODEL_PATH` | `OUTPUT_DIR / 'cellpose_GUV_model'` | Save path for the fine-tuned model |

#### 5-iii · Switch model

```python
USE_FINETUNED_MODEL = True   # True = use fine-tuned; False = revert to original model
```

Redefines `segment_guvs()` in-place so all downstream cells automatically use the new model.

#### 5-iv · Evaluate *(optional)*

Loads a held-out image/mask pair and runs both cyto3 and the fine-tuned model. Reports **Average Precision @ IoU = 0.5** and opens a napari comparison viewer with both predictions and the ground truth.

---

### Section 6 · Headless batch

Processes all `.vsi` files in `INPUT_DIR` without opening napari. Uses the current `segment_guvs()` definition (which reflects `CHOSEN_MODEL` and, if Section 5-iii was run, the fine-tuned model). Skips files where no GUVs are detected and logs errors without stopping the loop.

Outputs per-file `{stem}_measurements.csv` files and accumulates results into `all_dfs` for Section 7.

### Section 7 · Combined results

Merges results from Section 4c (interactive batch) and Section 6 (headless batch) — whichever was run — into `combined` and saves `all_GUV_measurements.csv`. Handles the case where only one of the two sections was run.

Produces `all_files_boxplots.png` — per-file boxplots with individual data points overlaid, for:
- GFP membrane mean intensity
- mCherry membrane mean intensity
- GFP membrane/lumen ratio
- mCherry membrane/lumen ratio

---

## Cellpose: model choice and parameter tuning

### cyto3 vs Cellpose-SAM (cpsam)

Cellpose offers two relevant pretrained models. The pipeline defaults to **cpsam**.

| | cyto3 | Cellpose-SAM (cpsam) |
|---|---|---|
| Architecture | U-Net (Cellpose native) | U-Net backbone initialised from a SAM encoder |
| Training data | ~700 k cell images across many cell types | ~1 M diverse images including cyto3 data |
| Size invariance | Weak — relies on `diameter` to normalise input | Strong — trained across 7.5–120 px; `diameter` is ignored |
| `diameter` parameter | Required for best results; pass `CP_DIAMETER` or `None` to auto-estimate | Not used; pass `None` or omit entirely |
| Cell probability threshold default | −2.0 (cyto3 convention) | 0.0 (cpsam convention) |
| Inference speed | Fast | ~2–3× slower than cyto3 |
| Typical accuracy on round objects | Good | Better on variable sizes and low-contrast boundaries |
| Download size | ~50 MB (bundled with Cellpose) | ~400 MB (downloaded on first use) |

**When to choose cyto3:**
- Speed is critical (large batches without GPU)
- GUVs are all similar in size and you can set `CP_DIAMETER` precisely
- cpsam produces too many false positives at `cellprob_threshold = 0`

**When to choose cpsam (default):**
- GUVs vary in size across images
- Low-contrast or dim GUVs that cyto3 misses
- No need to estimate or set `CP_DIAMETER`

### Parameter reference

#### `CP_DIAMETER` (cyto3 only)

The expected object diameter in pixels. Cellpose rescales the image so objects appear this size before running the neural network.

- `None`: auto-estimated from the image (slower; printed mean diameter is logged — copy it to `CP_DIAMETER` for future runs)
- Too small: large GUVs get clipped and may be split or missed
- Too large: small GUVs may be merged or missed
- Recommended: run once with `None`, note the printed mean, set `CP_DIAMETER` to that value

Not used by cpsam.

#### `CP_FLOW_THRESHOLD`

Controls how strictly Cellpose enforces gradient-flow consistency when assembling masks from the predicted flow field. It is a post-processing filter applied after the neural network output.

- Range: typically 0.0–1.0 (Cellpose default: 0.4)
- **Lower** → more permissive → keeps objects with noisier flows → more detections, more false positives
- **Higher** → stricter → discards objects with inconsistent flows → fewer detections, fewer false positives
- Pipeline default: `0` (maximally permissive) — appropriate because GUVs are well-defined rings and false positives can be removed in napari
- If you are getting many spurious detections in clean images, try raising to `0.3`–`0.4`

#### `CP_CELLPROB_THRESHOLD`

Threshold on the predicted cell-probability map. Pixels above this value are considered foreground before the flow-assembly step.

- **cyto3 convention**: typical range −6 to 0; default in the Cellpose GUI is −2
- **cpsam convention**: typical range −3 to 2; default is 0
- **Lower** → larger foreground region → detects dimmer and smaller GUVs, but increases background bleed-in
- **Higher** → smaller foreground region → only the brightest, most certain detections survive
- Pipeline default: `0` — a safe middle ground for cpsam. For cyto3 you will usually need `−2` or lower
- Troubleshooting: if GUVs are missed despite good signal, lower to `−2` or `−4`. If background patches are segmented, raise toward `+1` or `+2`

#### `CP_MIN_AREA_PX` / `MINIMUM_PHYSICAL_DIAMETER_UM`

`min_size` passed to Cellpose: any connected component smaller than this area (in pixels²) is discarded after mask assembly. Computed from the physical minimum diameter so it scales with pixel size.

- Increase `MINIMUM_PHYSICAL_DIAMETER_UM` to suppress small debris
- Decrease it to capture very small GUVs (but expect more noise objects)

#### `MEMBRANE_THICKNESS_PX`

Number of binary erosion steps applied inward from each GUV mask boundary to define the membrane ring. The lumen is the remaining eroded interior.

- Too thin: membrane ring only captures the outer edge; misses the inner leaflet; lumen bleeds into the ring
- Too thick: ring overlaps the lumen; lumen ROI becomes very small or disappears for small GUVs
- Verify using the ring-annotation figure in Section 4c: the orange overlay should match the bright GFP ring in the raw image
- Typical starting point: 5–10 px for images at 15 px/µm resolution
- GUVs too small to survive erosion (lumen disappears) are assigned the full mask as membrane with an empty lumen

### Typical tuning sequence

1. Run Section 4c on one representative file with `CP_DIAMETER = None`, `CP_FLOW_THRESHOLD = 0`, `CP_CELLPROB_THRESHOLD = 0`, `CHOSEN_MODEL = 'cpsam'`
2. Inspect the napari mask layer — note missed GUVs (→ lower `CP_CELLPROB_THRESHOLD`) or spurious detections (→ raise `CP_CELLPROB_THRESHOLD` or use napari to delete)
3. Check the segmentation preview figure — confirm mask IDs match expected GUV count
4. Check the ring-annotation figure — confirm the orange membrane ring matches the visible GFP ring; adjust `MEMBRANE_THICKNESS_PX` if needed
5. If switching to cyto3: copy the printed mean diameter to `CP_DIAMETER`, change `CP_CELLPROB_THRESHOLD` to `−2` as a starting point
6. Once satisfied, run Section 4c on the full `FILE_INDICES` list

---

## Output files

| File | Produced by | Description |
|---|---|---|
| `{stem}_segmentation_preview.png` | Section 4c | 4-panel: GFP \| mCherry \| masks with IDs \| composite overlay |
| `{stem}_ring_preview.png` | Section 4c | 4-panel membrane/lumen ring annotation |
| `{stem}_ratio_scatter.png` | Section 4c | GFP vs mCherry membrane/lumen ratio per GUV |
| `{stem}_lumen_ratio.png` | Section 4c | Background-subtracted lumen scatter coloured by GFP/mCherry ratio |
| `{stem}_measurements.csv` | Sections 4c, 6 | Per-GUV measurements for one file |
| `batch_interactive_measurements.csv` | Section 4c | All files from the interactive batch combined |
| `batch_interactive_lumen_ratio.png` | Section 4c | 3-panel combined plot: ratio scatter \| violin distributions \| ratio vs diameter |
| `batch_interactive_report.html` | Section 4c | Self-contained HTML report with all figures embedded as base64 |
| `training_data/{stem}_img.tif` | Section 5-i | GFP image saved for fine-tuning |
| `training_data/{stem}_masks.tif` | Section 5-i | Corrected mask saved for fine-tuning |
| `cellpose_GUV_model` | Section 5-ii | Fine-tuned Cellpose model weights |
| `all_GUV_measurements.csv` | Section 7 | All files from interactive and/or headless batch combined |
| `all_files_boxplots.png` | Section 7 | Per-file boxplots of key metrics |

---

## Measurement DataFrame columns

Each per-file CSV contains these columns. The interactive batch (4c) also adds `file_idx`. Both interactive (4c) and headless (6) batches produce the three background-subtracted lumen columns:

| Column | Unit | Description |
|---|---|---|
| `file` | — | VSI file stem |
| `file_idx` | — | Index into `vsi_files` *(interactive batch only)* |
| `label` | — | Integer GUV ID (1-based, assigned by Cellpose) |
| `centroid_x_px` | px | X centroid of the GUV mask |
| `centroid_y_px` | px | Y centroid |
| `area_px2` | px² | Total GUV area |
| `diameter_px` | px | Equivalent circular diameter: `2 × √(area/π)` |
| `GFP_membrane_mean` | raw counts | Mean GFP intensity in the membrane ring |
| `GFP_membrane_integrated` | raw counts | Summed GFP intensity in the membrane ring |
| `GFP_lumen_mean` | raw counts | Mean GFP intensity in the lumen |
| `GFP_lumen_integrated` | raw counts | Summed GFP intensity in the lumen |
| `GFP_mem_lumen_ratio` | — | `GFP_membrane_mean / (GFP_lumen_mean + ε)` |
| `mCherry_membrane_mean` | raw counts | Mean mCherry intensity in the membrane ring |
| `mCherry_membrane_integrated` | raw counts | Summed mCherry in the membrane ring |
| `mCherry_lumen_mean` | raw counts | Mean mCherry intensity in the lumen |
| `mCherry_lumen_integrated` | raw counts | Summed mCherry in the lumen |
| `mCherry_mem_lumen_ratio` | — | `mCherry_membrane_mean / (mCherry_lumen_mean + ε)` |
| `GFP_lumen_bgsub` | raw counts | `GFP_lumen_mean − GFP_bulk` clipped to 0 |
| `mCherry_lumen_bgsub` | raw counts | `mCherry_lumen_mean − mCherry_bulk` clipped to 0 |
| `GFP_mCherry_lumen_ratio` | — | `GFP_lumen_bgsub / (mCherry_lumen_bgsub + ε)` |

**Background (`_bulk`) estimation:** the mode of the histogram of all pixels outside every GUV mask (512 bins). This samples the external aqueous phase, which dominates the pixel-count distribution and represents the autofluorescence + detector offset baseline.

---

## Helper functions

### `read_vsi(vsi_path)`
Reads Series 0 of a CellSens VSI file using aicsimageio + BioFormats.  
Returns `(data, channel_names)` where `data` is a `(C, Y, X)` uint16 array.  
Raises `AssertionError` if Series 0 does not contain exactly 2 channels.

VSI series structure assumed:
- S=0: GFP + mCherry (2 planes, 2048×2048) — loaded
- S=1: Brightfield (1 plane) — ignored
- S=2: Macro image (3 planes, 512×512) — ignored

### `split_channels(data, channel_names, gfp_idx, mcherry_idx)`
Splits a CYX array by explicit channel index.  
Prints the full channel listing with assignments so you can verify every run.  
Returns `(gfp_img, mcherry_img)` as 2-D arrays.

### `normalise(img)`
Percentile stretch (1st–99th) to the 0–255 float32 range.  
Used for both Cellpose input and display overlays.

### `segment_guvs(gfp_img, diameter, flow_threshold, cellprob_threshold, model_type, use_gpu, min_size)`
Runs Cellpose on the normalised GFP image.  
Returns an `int32` label mask (0 = background).  
When `diameter=None`, prints the mean estimated diameter so you can set `CP_DIAMETER` for future runs.  
Redefined in Section 5-iii to use the fine-tuned model when `USE_FINETUNED_MODEL = True`.

### `make_membrane_lumen_masks(masks, membrane_thickness_px)`
Erodes each GUV binary mask inward by `membrane_thickness_px` steps using a 4-connected structuring element.  
Returns two dicts `{label → bool array}`: membrane ring (original minus eroded) and lumen (eroded interior).  
GUVs too small to survive full erosion retain the full mask as membrane with an empty lumen.

### `run_quantification(gfp_img, mcherry_img, masks, membrane_thickness_px)`
Calls `make_membrane_lumen_masks`, then calls `quantify()` for each channel.  
Returns `(df, membrane_masks, lumen_masks)`.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| No GUVs detected | `CP_CELLPROB_THRESHOLD` too high | Lower to `−2`, `−3`, or `−4` |
| Too many spurious detections | `CP_CELLPROB_THRESHOLD` too low | Raise toward `+1`; also check `CP_FLOW_THRESHOLD` (try `0.3`) |
| GUVs split into fragments | `CP_DIAMETER` too small (cyto3) | Increase `CP_DIAMETER` or use cpsam |
| Multiple GUVs merged into one | `CP_DIAMETER` too large (cyto3) | Decrease `CP_DIAMETER` |
| Small debris detected | `MINIMUM_PHYSICAL_DIAMETER_UM` too small | Increase to exclude debris size |
| Ring too thin / misses membrane | `MEMBRANE_THICKNESS_PX` too small | Increase by 2–3 px and re-check Section 4b |
| Ring bleeds into lumen | `MEMBRANE_THICKNESS_PX` too large | Decrease and re-check Section 4b |
| `GFP_mem_lumen_ratio` is NaN | GUV too small — lumen empty after erosion | Lower `MEMBRANE_THICKNESS_PX` or raise `MINIMUM_PHYSICAL_DIAMETER_UM` |
| Channels appear swapped | `GFP_CHANNEL_IDX` / `MCHERRY_CHANNEL_IDX` wrong | Swap the two values in Section 0b |
| napari window does not appear | Wrong kernel or missing Qt | Ensure `guv_env` kernel is selected; run `pip install pyqt5` |
| `AssertionError: Expected 2 channels` | File has unexpected series structure | Print `img.scenes` and adjust `S=` in `read_vsi` or the channel indices |

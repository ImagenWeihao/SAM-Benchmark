# SAM Benchmark -- MATLAB Demo

**SAM** (Smart Adaptive Microscopy) is a framework for evaluating how intelligently a microscopy system can allocate its imaging resources. This benchmark tests four levels of decision-making -- from naive random sampling up to LLM-guided adaptive control -- so their performance can be directly compared.

This repo currently implements **Levels 1, 2, and 3**, a standalone **Nucleus Normality Test** (Option 5), a same-zoom **Segmentation Validation** pipeline (Option 6), a cross-zoom **10X -> 40X Segmentation Validation** pipeline (Option 7) with automatic stitching, stage-based alignment, and NCC refinement, plus a batch-capable **GUI front-end** (`run_benchmark_gui`) with an optional real-time ND2 simulator. Level 4 (LLM-adaptive SAM) is in development.

---

## Entry Points

```
>> run_benchmark       % CLI -- numeric option menu
>> run_benchmark_gui   % GUI -- file/folder pickers, optional simulator
```

### CLI menu

```
[1]  Level 1 -- Random Sampling (Otsu segmentation)
[2]  Level 2 -- Circularity-Guided (rule-based)
[3]  Level 3 -- ML-Guided (UNet)
[4]  Run ALL three and compare
[5]  Nucleus Normality Test (405nm DAPI, MIP)
[6]  Segmentation Validation vs Ground Truth (PR curve + accuracy histogram)
[7]  10X -> 40X GT Validation (stitched, overlap crop, NCC-refined)
```

---

## Options 1-4 -- The Four Benchmark Levels

| Level | Selection Strategy | Segmentation | Status |
|-------|--------------------|--------------|--------|
| 1 | Random sampling | **Otsu + connected components** (non-if-then baseline) | Done |
| 2 | Circularity threshold | Adaptive threshold + morphology + watershed + circularity filter | Done |
| 3 | UNet confidence ranking | StarDist / Cellpose / MATLAB UNet | Done |
| 4 | LLM-adaptive | LLM reasons over image history | In development |

L1 and L2 now use **different** segmentation backbones: L1 is a pure Otsu baseline (no cleanup, no filter, no watershed) to establish a clean lower-tier comparison against L2's engineered pipeline. In earlier versions L1 and L2 shared `segment_nuclei.m` -- this is no longer the case.

All levels follow the same 5-step workflow for fair comparison:

1. Coarse scan -- detect all nuclei in the field of view
2. Subsample -- build a candidate pool (same budget across all levels)
3. Select targets -- each level uses its own strategy
4. Z-stack acquisition -- image selected cells across Z planes
5. Morphology characterisation -- measure shape, texture, and focus quality

Option **[4]** runs L1/L2/L3 sequentially on the same input and emits a side-by-side comparison.

### Shared parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_scan_cells` | 100 | Max cells in coarse scan pool |
| `n_zstack_cells` | 10 | Cells selected for Z-stack |
| `nucleus_channel` | 1 | Channel index for DAPI |
| `z_planes` | 7 | Z planes per cell (L1/L2) |
| `z_planes_highres` | 15 | Z planes per cell (L3) |
| `circularity_threshold` | 0.85 | L2 minimum circularity to qualify |
| `unet_model` | `'stardist'` | L3: `stardist` / `cellpose` / `nuclear_seg` / `matlab` / `mock` |
| `unet_threshold` | 0.75 | L3 minimum confidence score |

---

## Option 5 -- Nucleus Normality Test

Two-stage pipeline for nucleus morphology classification from the 405nm DAPI channel, operating on a **MIP** of the Z-stack rather than a single plane.

**Stage 1 -- Segmentation:** `segment_nuclei_normality` runs adaptive threshold + watershed with a second pass that catches small micronuclei.

**Stage 2 -- Classification** (set via `params.normality_mode`):

| Class | Morphologies | Signal |
|-------|-------------|--------|
| `normal` | Round, compact, single (A) | High prob, low shape variance |
| `abnormal_shape` | Multi-lobular, blebbing (B, F) | Low circularity or low solidity |
| `abnormal_count` | Binucleated, polyploid, micronuclei (C, D, E) | Area outlier |

### MIP parameters

```matlab
params.normality_use_mip         = true;
params.normality_series          = 1;     % first series of target FOV
params.normality_series_end      = 15;    % last series (= number of Z planes)
params.normality_nucleus_channel = 2;     % 405nm index in this .nd2 channel order
params.normality_mode            = 'stardist';
params.normality_circ_threshold  = 0.80;
params.normality_area_ratio_hi   = 2.5;
params.normality_area_ratio_lo   = 0.30;
```

---

## Option 6 -- Same-Zoom Segmentation Validation

Quantitative evaluation of each segmentation method against a manually annotated ground-truth mask and centroid table exported from NIS-Elements.

### Inputs

- **GT mask** -- multi-frame TIFF (Frame 0 = nucleus label mask, uint16)
- **GT annotations (optional)** -- CSV with columns `cell_id, morphology, nucleus_state, centroid_x, centroid_y, area`. Omit to disable ROI cropping; scoring then runs on the full image.
- **nd2 file** -- raw acquisition; correct series + channel must match the mask.

### Key parameters

```matlab
params.gt_mask_path       = 'path\to\masks_FOV1.tif';
params.gt_csv_path        = 'path\to\annotations_FOV1.csv';  % optional
params.gt_nd2_path        = 'path\to\raw.nd2';
params.gt_series          = 1;     % start series of target FOV
params.gt_series_end      = 15;    % end series (for MIP)
params.gt_use_mip         = true;  % MIP across Z slices before segmenting
params.gt_nucleus_channel = 2;
params.gt_roi_enabled     = true;  % interactive ROI draw on the 40X image
params.gt_export_centroids_csv = true;
```

### Methods compared

| Method | Segmentation backbone |
|--------|----------------------|
| Otsu (L1) | `graythresh` + `bwlabel` -- pure textbook baseline, no cleanup |
| Rule-based (L2) | Adaptive threshold + morphology + watershed + circularity filter (>= 0.60) |
| StarDist (L3) | Pretrained `2D_versatile_fluo` via Python bridge |

### Metrics

- **AP** (Average Precision): area under the PR curve swept over IoU thresholds 0.05-0.95 (step 0.05). Replaces the old pixel-area-score ROC AUC.
- **F1 @ IoU=0.5**: `2PR/(P+R)` from TPs matched at IoU >= 0.5 with ALL predictions (no score gating). Replaces the previous `f1_at_half` which used a misleading score threshold.
- **Precision-Recall curve**: Panel 4 sweeps IoU 0.05-0.95, highlights the IoU=0.5 operating point with a star marker.
- **Detection Rate (% of GT)**: Panel 5 shows GT = 100% baseline, Predicted/GT%, Matched/GT%.

### Main output figure (6 panels)

1-3. Per-method overlays -- predicted edges in method colour + GT edges in white, title shows AP / F1@IoU=0.5 / matched
4. Precision-Recall curve (IoU sweep 0.05-0.95) with IoU=0.5 operating point starred
5. Detection Rate (% of GT) bar chart
6. 40X + GT boundaries reference

### Accuracy-histogram figure (8 panels, when `gt_export_centroids_csv = true`)

1-3. Per-method detection-rate bars across IoU thresholds 0-1 (step 0.1) with IoU=0.5 marker + TP/FP/FN/Prec/Rec/F1/AP annotation box
4. All-method overlay line plot (detection rate vs IoU)
5-6. Paired accuracy bar at IoU=0.5 (L1/GT, L2/GT, L3/GT)
7-8. Grouped TP/FP/FN bars per method

### Centroid CSV export

When `gt_export_centroids_csv = true`, three files are written (one per method):

```
centroids_L1_<ts>.csv
centroids_L2_<ts>.csv
centroids_L3_<ts>.csv
```

Schema:

| Column | Description |
|--------|-------------|
| `cell_id` | Prediction label in the seg mask |
| `centroid_x_px`, `centroid_y_px` | Centroid in image pixels |
| `centroid_x_um`, `centroid_y_um` | Same in micrometres |
| `matched_gt_id` | GT cell ID matched at IoU >= 0.5 (0 = unmatched) |
| `iou_score` | Best IoU against any GT cell |

Rows with `matched_gt_id = 0` are false positives.

---

## Option 7 -- Cross-Zoom (10X -> 40X) Segmentation Validation

Segment a 10X overview mosaic and validate against a 40X ground-truth mask -- the two images are at different zooms and may come from different acquisition sessions.

### Workflow

1. **Stitch 10X mosaic** -- read `PlanePositionX/Y` from Bio-Formats OME metadata for every 10X series, invert per the Nikon `cameraTransformationMatrix` (stage +X -> image LEFT, stage +Y -> image UP), paste tiles onto a shared canvas in um coordinates.
2. **Load 40X FOVs + MIPs** -- Max Intensity Projection across the 15 Z planes per FOV; read each FOV's stage centre.
3. **Place 40X GT mask** in the 10X canvas -- downsample the native 40X label mask by `px_40x / px_10x ~= 0.248` using `'nearest'` (preserves integer labels).
4. **Crop overlap region** -- the area where both 10X image and 40X GT mask have content.
5. **Large-template NCC refinement** -- `normxcorr2` with the entire 10X segmentation mask as a single template against the padded GT mask, search radius `+/- params.align_ncc_search_um`. Multi-cell correspondences uniquely constrain the global shift.
6. **Optional ROI** -- when `gt_roi_enabled = true`, user draws a scoring rectangle on the NCC-refined 10X+GT overlay.
7. **Compute IoU-based PR curve + centroid CSVs** against the 10X segmentation.

### Key parameters

```matlab
params.nd2_10x_path               = 'path\to\10x_overview.nd2';
params.nd2_40x_path               = 'path\to\40x_highres.nd2';
params.channel_10x                = 1;
params.channel_40x                = 2;
params.n_fov_10x                  = 9;
params.n_fov_40x                  = 5;
params.z_planes_per_40x_fov       = 15;
params.align_ncc_search_um        = 100;   % NCC search radius
params.show_alignment_diagnostics = false; % true = pop up pre/post NCC diagnostic plots
params.gt_export_centroids_csv    = true;  % export per-method centroid CSVs
params.gt_roi_enabled             = true;  % interactive ROI on refined overlay
```

`show_alignment_diagnostics` is **off by default** to keep runs headless. Set true to see the pre-NCC / post-NCC / ROI-cropped overlap figures.

### Example result (WellA05 overlap, 20 GT nuclei)

| Method | Predicted | Matched @ IoU=0.5 | Recall | Precision | F1 | AP |
|--------|-----------|-------------------|--------|-----------|------|------|
| Otsu (L1) | 21 | 15/20 | 0.75 | 0.71 | 0.73 | 0.465 |
| Rule-based (L2) | 19 | 17/20 | 0.85 | 0.89 | 0.87 | 0.521 |
| StarDist (L3) | 21 | 19/20 | 0.95 | 0.90 | 0.93 | 0.475 |

L1 Otsu establishes a clear lower-tier baseline: more FPs and fewer matched GT than either L2 or L3. L2's circularity filter yields the highest precision. L3 achieves the highest recall with competitive precision. AP ordering can invert relative to F1 because PR curves from over-predicting methods sit lower overall.

### Output artefacts

- `overlap_diagnostic_initial_<ts>.png` -- 10X crop + GT mask at stage-based position (only when diagnostics enabled)
- `overlap_diagnostic_refined_<ts>.png` -- same after NCC refinement (only when diagnostics enabled)
- `validation_10x_40x_<ts>.png` -- 6-panel result: overlays + PR curve + Detection Rate bars + GT reference
- `accuracy_histogram_10x_<ts>.png` -- 8-panel detection-rate-vs-IoU summary (same layout as Option 6)
- `centroids_10x_L1_<ts>.csv`, `centroids_10x_L2_<ts>.csv`, `centroids_10x_L3_<ts>.csv`
- `results_validation_10x.mat` -- full result struct for downstream analysis

---

## GUI Front-End (`run_benchmark_gui`)

Three-panel configuration window with Start/Stop controls and a rolling log pane. Wraps all CLI options except interactive-only ones.

| Panel | Controls |
|-------|----------|
| 1. Input .nd2 | Single file / Folder radio, Browse, **Simulator mode** + interval (s) |
| 2. Ground Truth (optional) | Enable checkbox, GT mask .tif, GT CSV, GT .nd2, L1/L2/L3 level toggles |
| 3. Result Folder | Output folder (blank = auto `logs\<stem>_<timestamp>\`) |

### Simulator mode

When the input is a folder, files are processed one at a time with a user-set delay between them, simulating a real-time acquisition feed. Off by default; files are processed sequentially with no delay. Driven by a MATLAB `timer` object; Stop cleanly halts the timer and drains the queue.

### Example workflow -- GUI + GT Validation

1. Launch: `run_benchmark_gui`
2. Pick a single `.nd2` or a folder of `.nd2` files
3. Enable *GT benchmarking*, fill in:
   - `GT mask (.tif)` -- multi-frame TIFF with Frame 0 = label mask
   - `GT CSV` -- annotation table (optional; disables ROI cropping if blank)
   - `GT .nd2` -- raw nd2 used to generate the GT mask
4. Choose a *Result folder* (or leave blank to auto-save to `logs\<stem>_<timestamp>\`)
5. Click **Start**

Output artefacts per input file are the same as the CLI:

- `validation_10x_40x_<ts>.png` -- 6-panel overlay + PR curve + % of GT bars
- `accuracy_histogram_10x_<ts>.png` -- 8-panel detection-rate vs IoU summary
- `centroids_10x_L1/L2/L3_<ts>.csv` -- per-method centroid tables

Example centroid CSV row for an unmatched prediction (`matched_gt_id = 0`):

```
cell_id, centroid_x_px, centroid_y_px, centroid_x_um, centroid_y_um, matched_gt_id, iou_score
15,      355.63,        411.82,        232.92,        269.73,        0,             0.000
```

### When to use the GUI vs CLI

- **GUI**: batch runs, demos, when you want non-interactive scoring, or when you want to simulate a real-time feed.
- **CLI** (`run_benchmark`): interactive ROI selection, Option 7's full alignment pipeline with ROI draw, or quick single-file debugging.

---

## Requirements

### MATLAB
- R2022a or later
- Image Processing Toolbox
- Deep Learning Toolbox (Level 3 matlab mode only)

### Bio-Formats
- https://www.openmicroscopy.org/bio-formats/downloads/
- Unzip to e.g. `D:\Matlab Files\toolbox\bfmatlab\`

### Python (Level 3 + Normality Test + Validation)
- Python 3.9 standalone (not Anaconda, to avoid DLL conflicts)
- Configure in MATLAB: `pyenv('Version','C:\...\python.exe','ExecutionMode','OutOfProcess')`
- Install: `pip install "numpy<2" "tensorflow-cpu==2.10.0" "stardist==0.8.5"`

---

## Setup

1. Copy all files to one folder, e.g. `D:\Matlab Files\Codes\SAM\`
2. Edit the paths at the top of `run_benchmark.m` (CLI) or launch the GUI and use Browse buttons.
3. `startup.m`:
```matlab
pyenv('Version', 'C:\Users\...\Python39\python.exe', 'ExecutionMode', 'OutOfProcess');
addpath('D:\Matlab Files\toolbox\bfmatlab');
addpath('D:\Matlab Files\Codes\SAM');
jar_file = 'D:\Matlab Files\toolbox\bfmatlab\bioformats_package.jar';
if ~any(strcmp(javaclasspath('-dynamic'), jar_file)); javaaddpath(jar_file); end
```
4. Run: `run_benchmark` or `run_benchmark_gui`

> **Note on .m file encoding:** keep all .m files pure ASCII. MATLAB on Windows reads .m files as CP-1252 (not UTF-8), so non-ASCII characters (em dashes, box drawing, arrows, etc.) can silently break the parser and prevent the file from being registered -- symptom: `which my_file` returns 'not found' even when the file is present.

---

## File Overview

```
run_benchmark.m                      <- CLI entry point (numeric menu)
run_benchmark_gui.m                  <- GUI entry point (file/folder pickers, simulator)
load_nd2.m                           <- .nd2 reader with MIP support + multi-file cache
segment_nuclei.m                     <- L2/L3 segmentation (adaptive + watershed + solidity)
segment_nuclei_l1.m                  <- L1 segmentation (Otsu + connected components)
acquire_zstack.m                     <- Z-stack acquisition (shared)
characterize_morphology.m            <- Morphology features (shared)
display_results.m                    <- Figures and comparison plots
run_level1.m                         <- Level 1: random sampling + Otsu segmentation
run_level2.m                         <- Level 2: rule-based circularity filter
run_level3.m                         <- Level 3: ML-guided (UNet)
run_nucleus_normality_test.m        <- Option 5: normality classification (MIP-based)
run_segmentation_validation.m       <- Option 6: same-zoom PR-curve validation
run_segmentation_validation_10x.m   <- Option 7: 10X->40X stitching + NCC + PR curve
find_fov_location.m                 <- NCC diagnostic: locate single 40X FOV in 10X mosaic
classify_cells_unet.m               <- UNet dispatcher
unet_stardist.m / unet_cellpose.m /
unet_nuclear_seg.m / unet_mock.m /
unet_matlab.m                        <- Model adapters
stardist_bridge.py                   <- Python bridge (flat array transfer)
build_nucleus_classifier_unet.m      <- Lightweight UNet architecture
train_nucleus_classifier.m           <- Training scaffold
```

---

## Tested With

- MATLAB R2022a
- Python 3.9.13 (standalone, OutOfProcess)
- TensorFlow-CPU 2.10.0 + NumPy 1.24 + StarDist 0.8.5
- Nikon .nd2 -- 4-channel, WellA12: 1952x1952 px, 0.33 um/px
- Nikon .nd2 -- 4-channel, WellA05: 2048x2044 px, 0.1625 um/px, 5 FOVs
- Nikon .nd2 -- 3-channel + BF, WellA05: 10X overview, 9 FOVs

---

## Recent Changes

- **L1 now uses Otsu** -- previously shared `segment_nuclei.m` with L2. New `segment_nuclei_l1.m` (and a local copy inside each validation file) provides a textbook baseline. L2 keeps the adaptive pipeline + circularity filter.
- **ROC replaced with Precision-Recall over IoU sweep** -- the old pixel-area-score ROC was uninformative for object detection. Panel 4 is now a proper PR curve; AP replaces AUC; F1 is computed at IoU=0.5 with all predictions, not at a score threshold.
- **Detection Count bars replaced with % of GT** -- bars are now GT=100% baseline, Predicted/GT%, Matched/GT%, with per-bar percentage annotations.
- **Cross-zoom centroid CSV export** -- Option 7 now emits `centroids_10x_L1/L2/L3_<ts>.csv` matching Option 6's existing export.
- **Option 7 accuracy histogram** -- 8-panel detection-rate-vs-IoU summary (same layout as Option 6), written as `accuracy_histogram_10x_<ts>.png`.
- **Normality test uses MIP** -- Option 5 now runs on a Z-MIP instead of plane 1, with `params.normality_use_mip/series/series_end/nucleus_channel` controls.
- **`solidity` added to `segment_nuclei`** -- required by the normality mock classifier and the CSV export.
- **Alignment diagnostic toggle** -- `params.show_alignment_diagnostics` gates the four Option 7 overlap-diagnostic figures (default off).
- **GUI front-end + ND2 simulator** -- `run_benchmark_gui.m` provides file/folder pickers, a `timer`-driven drip simulator, and per-level toggles.
- **ASCII-only .m files** -- all non-ASCII characters stripped from the repo's .m files to avoid a Windows-specific MATLAB parser issue where UTF-8 files silently fail to register.

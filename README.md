# SAM Benchmark — MATLAB Demo

**SAM** (Smart Adaptive Microscopy) is a framework for evaluating how intelligently a microscopy system can allocate its imaging resources. This benchmark tests four levels of decision-making — from naive random sampling up to LLM-guided adaptive control — so their performance can be directly compared.

This repo currently implements **Levels 1, 2, and 3**, a standalone **Nucleus Normality Test**, a same-zoom **Segmentation Validation** pipeline (Option 6), and a **cross-zoom 10X→40X Segmentation Validation** pipeline (Option 7) with automatic stitching, stage-based alignment, and NCC refinement. Level 4 (LLM-adaptive SAM) is in development.

---

## The Four Levels

| Level | Name | Selection Strategy | Status |
|-------|------|--------------------|--------|
| 1 | Random Sampling | Randomly pick cells from the full FOV | Done |
| 2 | Rule-Based | Pick cells by circularity threshold | Done |
| 3 | ML-Guided | StarDist/Cellpose UNet classifier | Done |
| 4 | SAM — LLM Adaptive | LLM reasons over image history | In development |

All levels follow the same 5-step workflow for fair comparison:

1. Coarse scan — detect all nuclei in the field of view
2. Subsample — build a candidate pool (same budget across all levels)
3. Select targets — each level uses its own strategy
4. Z-stack acquisition — image selected cells across Z planes
5. Morphology characterisation — measure shape, texture, and focus quality

---

## Menu Options

```
[1]  Level 1 -- Random Sampling
[2]  Level 2 -- Circularity-Guided (rule-based)
[3]  Level 3 -- ML-Guided (UNet)
[4]  Run ALL three and compare
[5]  Nucleus Normality Test (405nm DAPI)
[6]  Segmentation Validation vs Ground Truth (ROC)
[7]  10X -> 40X GT Validation (stitched, overlap crop)
```

---

## Segmentation Validation (Option 6)

Quantitative evaluation of each segmentation method against manually annotated ground truth masks and centroid annotations exported from NIS-Elements.

### Inputs
- **GT mask** — multi-frame TIFF (Frame 0 = nucleus label mask, uint16)
- **GT annotations** — CSV with columns: `cell_id`, `morphology`, `nucleus_state`, `centroid_x`, `centroid_y`, `area`
- **nd2 file** — raw acquisition; correct series and channel must be confirmed

### Key parameters
```matlab
params.gt_mask_path       = 'path\to\masks_FOV1.tif';
params.gt_csv_path        = 'path\to\annotations_FOV1.csv';
params.gt_nd2_path        = 'path\to\raw.nd2';
params.gt_series          = 5;   % Bio-Formats series index for target FOV
params.gt_z_plane         = 1;   % Z plane within that series
params.gt_nucleus_channel = 2;   % channel index for DAPI (file-dependent)
```

> **Note on series selection:** Nikon nd2 files with multiple FOVs are stored as separate Bio-Formats series. The correct series must be confirmed by overlaying GT centroids on the loaded image. Use the explorer script to check all series before running validation.

### Methods compared

| Method | Description |
|--------|-------------|
| MATLAB (L1/L2) | Adaptive threshold + watershed segmentation |
| Circularity-filtered (L2) | Same, filtered by circularity >= threshold |
| StarDist (L3) | Pretrained `2D_versatile_fluo` via Python bridge |

### ROC metric
Per-nucleus IoU matching at thresholds 0 to 1. Primary metric is **F1 at IoU=0.5**. AUC is also reported but note that pixel-based FPR compresses the ROC curve for sparse FOVs — F1 is more interpretable.

### Validated results (WellA05, FOV1, 405nm DAPI, n=50 GT nuclei)

| Method | Detected | F1@IoU=0.5 | AUC |
|--------|----------|------------|-----|
| MATLAB (L1/L2) | 50/50 | 0.720 | 0.031 |
| Circ.-filtered (L2) | 47/50 | 0.701 | 0.027 |
| StarDist (L3) | 34/50 | 0.381 | 0.006 |

**Key findings:**
- MATLAB adaptive segmentation outperforms pretrained StarDist on this dataset. The nuclear morphology (ring-like DAPI staining) differs from StarDist's training distribution.
- Circularity filtering (L2) adds little value when the population is already highly circular (mean=0.862).
- StarDist struggles with blebbing nuclei due to its star-convex polygon constraint, which cannot represent concave boundaries.
- Both methods miss out-of-focus nuclei in the top-left region and touching nuclei in the dense center cluster.

**Known segmentation limitations:**
- Raising adaptive threshold sensitivity to catch diffuse nuclei causes severe over-segmentation on images with uneven background illumination. Background subtraction would be needed first.
- Lowering watershed `imhmin` depth to split touching nuclei increases false splits on elongated/blebbing nuclei. These parameters require per-dataset tuning.

### Output figure (6 panels)
1. MATLAB overlay — white=GT boundary, coloured=predicted, red X=missed cells with ID
2. Circularity-filtered overlay — same format
3. StarDist overlay — same format
4. GT annotation map — coloured by morphology class
5. ROC curves — all three methods with AUC and F1@0.5 operating points
6. Detection rate by morphology class (normal / blebbing / micronuclei / unsure)

---

## Cross-Zoom Segmentation Validation (Option 7)

Segment a 10X overview mosaic and validate against a 40X ground-truth mask — the two images are at different zooms and may come from different acquisition sessions.

### Workflow
1. **Stitch 10X mosaic** — read `PlanePositionX/Y` from Bio-Formats OME metadata for every 10X series, invert per the Nikon `cameraTransformationMatrix` (stage +X → image LEFT, stage +Y → image UP), and paste tiles onto a shared canvas in µm coordinates.
2. **Load 40X FOVs + MIPs** — compute a Max Intensity Projection across the 15 Z planes per FOV, read each FOV's stage centre.
3. **Place 40X GT mask** in the 10X canvas — downsample the native 40X label mask by `px_40x / px_10x ≈ 0.248` using `'nearest'` (preserves integer labels), positioned at the physical stage location.
4. **Crop overlap region** — the area where both 10X image and 40X GT mask have content.
5. **Large-template NCC refinement** — stage coordinates drift between sessions (typically 30–50 µm). Run `normxcorr2` with the **entire 10X segmentation mask** as a single template against the padded GT mask, search ±100 µm. The multi-cell cluster provides many simultaneous correspondences which constrain the global shift uniquely; this is more robust than sliding small windows over single cells.
6. **Apply the corrected shift** (sign-inverted relative to NCC peak location) to the GT mask, then compute IoU-based ROC against the 10X segmentation.

### Key parameters
```matlab
params.nd2_10x_path          = 'path\to\10x_overview.nd2';
params.nd2_40x_path          = 'path\to\40x_highres.nd2';
params.channel_10x           = 1;    % DAPI channel in 10X file
params.channel_40x           = 2;    % DAPI channel in 40X file
params.n_fov_10x             = 9;    % 3x3 10X mosaic
params.n_fov_40x             = 5;    % five 40X FOVs
params.z_planes_per_40x_fov  = 15;
params.align_ncc_search_um   = 100;  % NCC search radius
```

### Validated results (WellA05, FOV1 overlap, 405nm DAPI)

| Method | Detected | GT in overlap | F1 @ IoU=0.5 | AUC |
|--------|----------|---------------|--------------|-----|
| MATLAB (L1/L2) | 28 | 22 | 0.000 | 0.100 |
| Circularity-filtered (L2) | 27 | 22 | **0.600** | 0.107 |
| StarDist (L3) | 29 | 22 | **0.773** | 0.143 |

**NCC refinement** found a global shift of (−7, −65) px = (−4.6, −42.6) µm between the stage-predicted position and the true physical location, with peak NCC = 0.511. Before the refinement, all three methods scored 0.0 on F1 — the entire validation depends on this correction step.

### Output figures
- `overlap_diagnostic_initial_<ts>.png` — 10X crop + GT mask at stage-based position
- `overlap_diagnostic_refined_<ts>.png` — same after NCC refinement (GT now overlays 10X nuclei)
- `validation_10x_40x_<ts>.png` — 6-panel result: three segmentation overlays + ROC + detection counts + GT reference

---

## Nucleus Normality Test (Option 5)

Two-stage pipeline for nucleus morphology classification from the 405nm DAPI channel.

**Stage 1 — Segmentation:** MATLAB adaptive threshold + watershed with a second micronuclei detection pass.

**Stage 2 — Classification** (set via `params.normality_mode`):

| Class | Morphologies | Signal |
|-------|-------------|--------|
| `normal` | Round, compact, single (A) | High prob, low shape variance |
| `abnormal_shape` | Multi-lobular, blebbing (B, F) | High radial distance CoV |
| `abnormal_count` | Binucleated, polyploid, micronuclei (C, D, E) | n_objects > 1 or area outlier |

**Validated on WellA12:** 15 nuclei, 67% normal, 27% abnormal shape, 7% abnormal count.

---

## Requirements

### MATLAB
- R2022a or later
- Image Processing Toolbox
- Deep Learning Toolbox (Level 3 matlab mode only)

### Bio-Formats
- Download: https://www.openmicroscopy.org/bio-formats/downloads/
- Unzip to e.g. `D:\Matlab Files\toolbox\bfmatlab\`

### Python (Level 3 + Normality Test + Validation)
- Python 3.9 from https://python.org/downloads/release/python-3913/
- Configure in MATLAB: `pyenv('Version', 'C:\path\to\python39.exe', 'ExecutionMode', 'OutOfProcess')`
- Install: `pip install "numpy<2" "tensorflow-cpu==2.10.0" "stardist==0.8.5"`

> R2022a supports Python 3.9 only. Use standalone Python (not Anaconda) with OutOfProcess mode to avoid DLL conflicts.

---

## Setup

**1.** Copy all files to one folder, e.g. `D:\Matlab Files\Codes\SAM\`

**2.** Edit paths at the top of `run_benchmark.m`:
```matlab
bf_path         = 'D:\Matlab Files\toolbox\bfmatlab';
params.nd2_path = 'C:\path\to\your\file.nd2';
```

**3.** Add to `startup.m`:
```matlab
pyenv('Version', 'C:\Users\...\Python39\python.exe', 'ExecutionMode', 'OutOfProcess');
addpath('D:\Matlab Files\toolbox\bfmatlab');
addpath('D:\Matlab Files\Codes\SAM');
jar_file = 'D:\Matlab Files\toolbox\bfmatlab\bioformats_package.jar';
if ~any(strcmp(javaclasspath('-dynamic'), jar_file)); javaaddpath(jar_file); end
```

**4.** Run: `run_benchmark`

---

## Parameters

### Shared
| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_scan_cells` | 100 | Max cells in coarse scan pool |
| `n_zstack_cells` | 10 | Cells selected for Z-stack |
| `nucleus_channel` | 1 | Channel index for DAPI |
| `z_planes` | 7 | Z planes per cell (L1/L2) |
| `z_planes_highres` | 15 | Z planes per cell (L3) |

### Level 2
| Parameter | Default | Description |
|-----------|---------|-------------|
| `circularity_threshold` | 0.85 | Minimum circularity to qualify |

### Level 3
| Parameter | Default | Description |
|-----------|---------|-------------|
| `unet_model` | `'mock'` | `stardist` / `cellpose` / `nuclear_seg` / `matlab` / `mock` |
| `unet_threshold` | 0.75 | Minimum confidence score |

### Normality Test
| Parameter | Default | Description |
|-----------|---------|-------------|
| `normality_mode` | `'mock'` | `stardist` / `mock` |
| `normality_circ_threshold` | 0.80 | Below = abnormal_shape |
| `normality_area_ratio_hi` | 2.5 | Above = polyploid |
| `normality_area_ratio_lo` | 0.30 | Below = micronucleus |

---

## File Overview

```
run_benchmark.m                ← Entry point
load_nd2.m                     ← .nd2 reader (Bio-Formats, multi-file cache)
segment_nuclei.m               ← Nucleus detection (shared, pixel-size aware)
acquire_zstack.m               ← Z-stack acquisition (shared)
characterize_morphology.m      ← Morphology features (shared)
display_results.m              ← Figures and comparison plots
run_level1.m                   ← Level 1: random sampling
run_level2.m                   ← Level 2: circularity-guided
run_level3.m                   ← Level 3: ML-guided (UNet)
run_nucleus_normality_test.m   ← Nucleus normality classification
run_segmentation_validation.m  ← Ground truth ROC validation (Option 6, same-zoom)
run_segmentation_validation_10x.m ← Cross-zoom 10X→40X validation (Option 7, stitching + NCC refinement)
find_fov_location.m            ← NCC-based location of a single 40X FOV in the 10X mosaic (diagnostic)
export_overlap_tif.py          ← Python reference: stage + cameraTransformationMatrix alignment → overlap TIFs
classify_cells_unet.m          ← UNet model dispatcher
unet_stardist.m                ← StarDist adapter
unet_cellpose.m                ← Cellpose adapter
unet_nuclear_seg.m             ← nuclearSegmentator adapter
unet_mock.m                    ← Rule-based mock
unet_matlab.m                  ← MATLAB Deep Learning Toolbox adapter
stardist_bridge.py             ← Python bridge (flat array transfer)
build_nucleus_classifier_unet.m← Lightweight UNet architecture
train_nucleus_classifier.m     ← Training scaffold
```

---

## Tested With

- MATLAB R2022a
- Python 3.9.13 (standalone, OutOfProcess)
- TensorFlow-CPU 2.10.0 + NumPy 1.24 + StarDist 0.8.5
- Nikon .nd2 — 4-channel, WellA12: 1952x1952px, 0.33µm/px
- Nikon .nd2 — 4-channel, WellA05: 2048x2044px, 0.1625µm/px, 5 FOVs
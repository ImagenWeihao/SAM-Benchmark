# SAM Benchmark — MATLAB Demo

**SAM** (Sequential Adaptive Microscopy) is a framework for evaluating how intelligently a microscopy system can allocate its imaging resources. This benchmark tests four levels of decision-making — from naive random sampling up to LLM-guided adaptive control — so their performance can be directly compared.

This repo currently implements **Levels 1, 2, and 3**, plus a standalone **Nucleus Normality Test** using StarDist. Level 4 (LLM-adaptive SAM) is in development.

---

## The Four Levels

| Level | Name | Selection Strategy | Status |
|-------|------|--------------------|--------|
| 1 | Random Sampling | Randomly pick cells from the full FOV | ✅ Done |
| 2 | Rule-Based | Pick cells by circularity threshold (if-then) | ✅ Done |
| 3 | ML-Guided | UNet classifier scores candidates | ✅ Done |
| 4 | SAM — LLM Adaptive | LLM reasons over image history | 🔧 In development |

All levels follow the same 5-step workflow for fair comparison:

1. Coarse scan — detect all nuclei in the field of view
2. Subsample — build a candidate pool (same budget across all levels)
3. Select targets — each level uses its own strategy here
4. Z-stack acquisition — image selected cells across Z planes
5. Morphology characterisation — measure shape, texture, and focus quality

---

## Nucleus Normality Test (Option 5)

A standalone two-stage pipeline for nucleus morphology classification from the **405nm DAPI channel only**:

**Stage 1 — Segmentation:** MATLAB adaptive threshold + watershed detects all nuclei, with a second micronuclei pass for small satellite objects.

**Stage 2 — Classification:** StarDist `2D_versatile_fluo` pretrained model classifies each nucleus into one of three classes:

| Class | Morphologies | Key Signal |
|-------|-------------|------------|
| `normal` | Round, compact, single nucleus (A) | High prob, low shape_var |
| `abnormal_shape` | Multi-lobular, blebbing (B, F) | High shape_var from radial profile |
| `abnormal_count` | Binucleated, polyploid, micronuclei (C, D, E) | n_objects > 1 or area outlier |

Reference morphology diagram:
- A — Normal, B — Multi-lobular, C — Binucleated, D — Polyploid
- E — Micronuclei, F — Blebbing nucleus, G — Advanced substructure

**Validated result** on `WellA12_Channel405,488,561,640`:
- 15 nuclei detected (12 standard + 3 micronuclei)
- 10 normal (66.7%), 4 abnormal shape (26.7%), 1 abnormal count (6.7%)
- StarDist runtime: ~1.6s for 15 cells

---

## Requirements

### MATLAB
- **R2022a** or later
- Image Processing Toolbox
- Deep Learning Toolbox *(Level 3 matlab mode only)*

### Bio-Formats
Required to read `.nd2` files:
- Download: https://www.openmicroscopy.org/bio-formats/downloads/
- Unzip to e.g. `D:\Matlab Files\toolbox\bfmatlab\`

### Python (Level 3 + Normality Test)
- **Python 3.9** — install from https://python.org/downloads/release/python-3913/
- Configure in MATLAB once: `pyenv('Version', 'C:\path\to\python39.exe', 'ExecutionMode', 'OutOfProcess')`

Install StarDist dependencies (Command Prompt):
```bash
pip install "numpy<2" "tensorflow-cpu==2.10.0" "stardist==0.8.5"
```

> **Note for R2022a users:** Python 3.10+ is not supported. Use Python 3.9 standalone (not Anaconda) with `ExecutionMode OutOfProcess` to avoid DLL conflicts.

---

## Setup

**1.** Copy all files into one folder, e.g. `D:\Matlab Files\Codes\SAM\`

**2.** Edit the paths at the top of `run_benchmark.m`:
```matlab
bf_path         = 'D:\Matlab Files\toolbox\bfmatlab';
params.nd2_path = 'C:\path\to\your\file.nd2';
```

**3.** Add to `startup.m` for automatic setup every session:
```matlab
pyenv('Version', 'C:\Users\...\Python39\python.exe', 'ExecutionMode', 'OutOfProcess');
addpath('D:\Matlab Files\toolbox\bfmatlab');
addpath('D:\Matlab Files\Codes\SAM');
jar_file = 'D:\Matlab Files\toolbox\bfmatlab\bioformats_package.jar';
if ~any(strcmp(javaclasspath('-dynamic'), jar_file)); javaaddpath(jar_file); end
```

**4.** Run:
```matlab
run_benchmark
```

---

## Menu Options

```
[1]  Level 1 -- Random Sampling
[2]  Level 2 -- Circularity-Guided (rule-based)
[3]  Level 3 -- ML-Guided (UNet)
[4]  Run ALL three and compare
[5]  Nucleus Normality Test (405nm DAPI)
```

---

## Parameters

All parameters are in the **User Parameters** section of `run_benchmark.m`:

### Shared
| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_scan_cells` | 100 | Max cells in coarse scan pool |
| `n_zstack_cells` | 10 | Cells selected for Z-stack |
| `nucleus_channel` | 1 | Channel index for DAPI / 405nm |
| `z_planes` | 7 | Z planes per cell (L1/L2) |
| `z_planes_highres` | 15 | Z planes per cell (L3) |

### Level 2
| Parameter | Default | Description |
|-----------|---------|-------------|
| `circularity_threshold` | 0.85 | Minimum circularity to qualify |

### Level 3 — UNet
| Parameter | Default | Description |
|-----------|---------|-------------|
| `unet_model` | `'mock'` | Model: `stardist` / `cellpose` / `nuclear_seg` / `matlab` / `mock` |
| `unet_target` | `'round_nucleus'` | Target morphology description |
| `unet_threshold` | 0.75 | Minimum confidence score |

### Nucleus Normality Test
| Parameter | Default | Description |
|-----------|---------|-------------|
| `normality_mode` | `'mock'` | Mode: `stardist` / `mock` |
| `normality_circ_threshold` | 0.80 | Below = abnormal_shape candidate |
| `normality_solid_threshold` | 0.85 | Below = binucleation candidate |
| `normality_area_ratio_hi` | 2.5 | Above = polyploid candidate |
| `normality_area_ratio_lo` | 0.30 | Below = micronucleus candidate |

---

## UNet Model Options (Level 3 + Normality Test)

| Model | `params.unet_model` | Install | Notes |
|-------|-------------------|---------|-------|
| StarDist | `'stardist'` | `pip install stardist tensorflow-cpu==2.10.0` | **Recommended** for DAPI nuclei |
| Cellpose | `'cellpose'` | `pip install cellpose` | Good general segmentation |
| nuclearSegmentator | `'nuclear_seg'` | `pip install nuclear-segmentator` | Direct normality classification |
| MATLAB UNet | `'matlab'` | Deep Learning Toolbox | Needs trained `.mat` model |
| Mock | `'mock'` | None | Rule-based proxy, no install needed |

---

## Output

Each run creates a timestamped log folder next to the `.m` files:

```
SAM/logs/20260414_161952/
├── run_log.txt                    ← full console output
├── nucleus_raw.png                ← raw 405nm channel
├── Level1_<timestamp>.png         ← L1 results figure
├── Level2_<timestamp>.png         ← L2 results figure
├── Level3_<timestamp>.png         ← L3 results figure
├── Comparison_<timestamp>.png     ← side-by-side comparison
├── normality_test_<timestamp>.png ← normality classification figure
├── results_Level1.mat             ← full L1 results struct
├── results_Level2.mat             ← full L2 results struct
├── results_Level3.mat             ← full L3 results struct
├── results_normality.mat          ← normality test results
├── morphology_Level1.csv          ← per-cell morphology table
├── morphology_Level2.csv
├── morphology_Level3.csv          ← includes UNet score column
└── normality_classification.csv   ← per-cell class, confidence, reason
```

---

## File Overview

```
run_benchmark.m                ← Entry point — run this
load_nd2.m                     ← .nd2 file reader (Bio-Formats)
segment_nuclei.m               ← Nucleus detection (shared)
acquire_zstack.m               ← Z-stack acquisition (shared)
characterize_morphology.m      ← Morphology features (shared)
display_results.m              ← Figures and comparison plots
run_level1.m                   ← Level 1: random sampling
run_level2.m                   ← Level 2: circularity-guided
run_level3.m                   ← Level 3: ML-guided (UNet)
run_nucleus_normality_test.m   ← Nucleus normality classification
classify_cells_unet.m          ← UNet model dispatcher
unet_stardist.m                ← StarDist adapter (recommended)
unet_cellpose.m                ← Cellpose adapter
unet_nuclear_seg.m             ← nuclearSegmentator adapter
unet_matlab.m                  ← MATLAB Deep Learning Toolbox adapter
unet_mock.m                    ← Rule-based mock (no install needed)
stardist_bridge.py             ← Python bridge for StarDist
build_nucleus_classifier_unet.m← Lightweight UNet architecture builder
train_nucleus_classifier.m     ← Training scaffold for custom UNet
```

---

## Tested With

- MATLAB R2022a
- Python 3.9.13 (standalone, OutOfProcess mode)
- TensorFlow-CPU 2.10.0 + NumPy 1.24 + StarDist 0.8.5
- Nikon `.nd2` — 4-channel well-plate, 1952×1952px, 0.33µm/px
- File: `WellA12_Channel405,488,561,640_Seq0011.nd2`
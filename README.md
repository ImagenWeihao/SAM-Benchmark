# SAM Benchmark — MATLAB Demo

**SAM** (Sequential Adaptive Microscopy) is a framework for evaluating how intelligently a microscopy system can allocate its imaging resources. This benchmark tests four levels of decision-making, from naive random sampling up to LLM-guided adaptive control, so their performance can be directly compared.

This repo currently implements **Levels 1 and 2** (traditional baselines). Levels 3 and 4 are in development.

---

## The Four Levels

| Level | Name | Selection Strategy |
|-------|------|--------------------|
| 1 | Random Sampling | Randomly pick cells from the full FOV |
| 2 | Rule-Based | Pick cells by circularity threshold (if-then rule) |
| 3 | ML-Guided *(coming)* | Use a pretrained model to predict imaging targets |
| 4 | SAM — LLM Adaptive *(coming)* | LLM reasons over image history to guide acquisition |

All levels follow the same 5-step workflow so results are directly comparable:

1. Coarse scan — detect all nuclei in the field of view
2. Subsample — build a candidate pool (same budget across all levels)
3. Select targets — each level uses its own strategy here
4. Z-stack acquisition — image selected cells across Z planes
5. Morphology characterisation — measure shape, texture, and focus quality

---

## Requirements

- MATLAB with **Image Processing Toolbox**
- **Bio-Formats** MATLAB toolbox (`bfopen`) — for reading `.nd2` files
  - Download: https://www.openmicroscopy.org/bio-formats/downloads/
  - Unzip to a local folder, e.g. `D:\Matlab Files\toolbox\bfmatlab\`

---

## Setup

1. Copy all `.m` files into a single folder, e.g. `D:\Matlab Files\Codes\SAM\`
2. Open MATLAB and navigate (`cd`) to that folder
3. Edit the two paths at the top of `run_benchmark.m`:

```matlab
bf_path      = 'D:\Matlab Files\toolbox\bfmatlab';   % Bio-Formats location
params.nd2_path = 'C:\path\to\your\file.nd2';         % Your .nd2 image file
```

4. Run:
```matlab
run_benchmark
```

Bio-Formats is loaded automatically — no manual `addpath` needed.

---

## Parameters

All tunable parameters are in the **User Parameters** section of `run_benchmark.m`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_scan_cells` | 100 | Max cells in the coarse scan pool |
| `n_zstack_cells` | 10 | Number of cells selected for Z-stack |
| `circularity_threshold` | 0.85 | Level 2: minimum circularity to qualify (0–1) |
| `nucleus_channel` | 1 | Channel index for DAPI / 405 nm |
| `z_planes` | 7 | Z planes per cell (simulated if file is 2D-only) |
| `z_step_um` | 0.5 | Microns between Z planes |

---

## Output

Each run automatically creates a timestamped log folder next to the `.m` files:

```
SAM/
└── logs/
    └── 20260410_143022/
        ├── run_log.txt              ← full console output
        ├── nucleus_raw.png          ← raw nucleus channel image
        ├── results_Level1.mat       ← full results struct (MATLAB)
        ├── morphology_Level1.csv    ← per-cell morphology table
        └── Level1_20260410_143022.png  ← figure screenshot
```

The morphology CSV contains: Cell ID, area, circularity, aspect ratio, solidity, texture contrast, and focus quality score.

---

## File Overview

```
run_benchmark.m           ← Entry point — run this
load_nd2.m                ← Reads .nd2 files via Bio-Formats
segment_nuclei.m          ← Nucleus detection and measurement (shared)
acquire_zstack.m          ← Z-stack extraction or simulation (shared)
characterize_morphology.m ← Morphology feature extraction (shared)
run_level1.m              ← Level 1: random sampling
run_level2.m              ← Level 2: circularity-guided selection
display_results.m         ← Figures and comparison plots
```

---

## Tested With

- MATLAB R2022b or later
- Nikon `.nd2` files (multi-channel, well-plate format)
- Example file: `WellA12_Channel405,488,561,640` — 4 channels, 1952×1952 px, 0.33 µm/px

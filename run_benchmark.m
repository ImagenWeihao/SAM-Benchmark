%% SAM Benchmark Demo -- Level 1 & 2
% Entry point: loads .nd2 file and runs selected benchmark level
%
% Requirements:
%   - Bio-Formats MATLAB toolbox (bfopen): https://www.openmicroscopy.org/bio-formats/
%   - Image Processing Toolbox
%
% Usage:
%   >> run_benchmark
%
% Author: SAM Project
% Date:   2026

clear; clc; close all;

%% ── Bio-Formats Path Setup ───────────────────────────────────────────────
bf_path = 'D:\Matlab Files\toolbox\bfmatlab';

if ~isfolder(bf_path)
    error('SAM:pathNotFound', ...
          'bfmatlab folder not found at:\n  %s\nUpdate bf_path in run_benchmark.m', bf_path);
end
addpath(bf_path);

jar_file = fullfile(bf_path, 'bioformats_package.jar');
if ~exist(jar_file, 'file')
    error('SAM:jarNotFound', 'bioformats_package.jar not found in:\n  %s', bf_path);
end
existing_cp = javaclasspath('-dynamic');
if ischar(existing_cp); existing_cp = {existing_cp}; end
if ~any(strcmp(existing_cp, jar_file))
    javaaddpath(jar_file);
    fprintf('[INFO] Loaded: bioformats_package.jar\n');
end
if ~exist('bfopen', 'file')
    error('SAM:bfNotFound', 'bfopen.m not found after addpath. Check:\n  %s', bf_path);
end
fprintf('[INFO] Bio-Formats ready.\n');

%% ── User Parameters ──────────────────────────────────────────────────────
params.nd2_path        = 'C:\Users\chenw\Documents\Downloads\20251024_151823_127__WellA12_Channel405,488,561,640_Seq0011.nd2';

% Level 1 params
params.n_scan_cells    = 100;   % max cells in coarse scan pool
params.n_zstack_cells  = 10;    % cells selected for Z-stack

% Level 2 params
params.circularity_threshold = 0.85;  % 0-1; higher = more circular

% Channel indices (1-based: 405, 488, 561, 640 nm)
params.nucleus_channel = 1;
params.gfp_channel     = 2;
params.rfp_channel     = 3;
params.far_red_channel = 4;

% Z-stack params (used if real Z-stack absent in file)
params.z_planes   = 7;
params.z_step_um  = 0.5;

%% ── Logging Setup ────────────────────────────────────────────────────────
% Log folder is created next to this script file
script_dir          = fileparts(mfilename('fullpath'));
params.run_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
params.log_dir       = fullfile(script_dir, 'logs', params.run_timestamp);

if ~exist(params.log_dir, 'dir')
    mkdir(params.log_dir);
end
fprintf('[LOG]  Run log folder: %s\n', params.log_dir);

% Open log file
log_file = fullfile(params.log_dir, 'run_log.txt');
diary(log_file);
fprintf('[LOG]  Console output mirrored to: run_log.txt\n\n');

%% ── Select Level ─────────────────────────────────────────────────────────
fprintf('\n+----------------------------------------------+\n');
fprintf('|       SAM Benchmark Demo  --  Level 1 & 2   |\n');
fprintf('+----------------------------------------------+\n\n');
fprintf('Select benchmark level:\n');
fprintf('  [1]  Level 1 -- Random Sampling Baseline\n');
fprintf('  [2]  Level 2 -- Circularity-Guided Selection\n');
fprintf('  [3]  Run BOTH and compare side-by-side\n\n');

level_choice = input('Enter choice (1/2/3): ');

%% ── Load .nd2 File ───────────────────────────────────────────────────────
fprintf('\n[INFO] Loading .nd2 file...\n');
fprintf('       %s\n\n', params.nd2_path);

img_data = load_nd2(params.nd2_path, params);

% Save a copy of the nucleus channel image to log
nucleus_fig = figure('Visible', 'off');
imagesc(img_data.nucleus); colormap gray; axis image off;
title(sprintf('Nucleus channel (405nm) -- %s', params.run_timestamp), 'Interpreter', 'none');
exportgraphics(nucleus_fig, fullfile(params.log_dir, 'nucleus_raw.png'), 'Resolution', 150);
close(nucleus_fig);
fprintf('[LOG]  Raw nucleus image saved.\n');

%% ── Run Selected Level(s) ────────────────────────────────────────────────
switch level_choice
    case 1
        results_L1 = run_level1(img_data, params);
        display_results(results_L1, [], params);
        save_results(results_L1, [], params);

    case 2
        results_L2 = run_level2(img_data, params);
        display_results([], results_L2, params);
        save_results([], results_L2, params);

    case 3
        results_L1 = run_level1(img_data, params);
        results_L2 = run_level2(img_data, params);
        display_results(results_L1, results_L2, params);
        compare_levels(results_L1, results_L2, params);
        save_results(results_L1, results_L2, params);

    otherwise
        error('Invalid choice. Please enter 1, 2, or 3.');
end

fprintf('\n[DONE] Benchmark complete.\n');
fprintf('[LOG]  All outputs saved to: %s\n', params.log_dir);
diary off;


%% ── Helper: save results struct and CSV ─────────────────────────────────
function save_results(R1, R2, params)
% Save .mat results and morphology CSV for each level that was run.

levels = {};
if ~isempty(R1); levels{end+1} = R1; end
if ~isempty(R2); levels{end+1} = R2; end

for i = 1:numel(levels)
    R = levels{i};

    % Save full results struct
    mat_name = fullfile(params.log_dir, sprintf('results_Level%d.mat', R.level));
    save(mat_name, 'R');
    fprintf('[LOG]  Results saved: results_Level%d.mat\n', R.level);

    % Save morphology table as CSV
    csv_name = fullfile(params.log_dir, sprintf('morphology_Level%d.csv', R.level));
    fid = fopen(csv_name, 'w');
    fprintf(fid, 'CellID,Area_px2,Circularity,AspectRatio,Solidity,TextureContrast,FocusQuality\n');
    for k = 1:numel(R.morph)
        m = R.morph(k);
        fprintf(fid, '%d,%.2f,%.4f,%.4f,%.4f,%.6f,%.6f\n', ...
                m.cell_id, m.nucleus_area, m.circularity, ...
                m.aspect_ratio, m.solidity, m.texture_contrast, m.focus_quality);
    end
    fclose(fid);
    fprintf('[LOG]  Morphology CSV saved: morphology_Level%d.csv\n', R.level);
end
end

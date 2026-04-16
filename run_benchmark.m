%% SAM Benchmark Demo -- Levels 1, 2 & 3
% Entry point: loads .nd2 file and runs selected benchmark level(s)
%
% Requirements:
%   - Bio-Formats MATLAB toolbox (bfopen): https://www.openmicroscopy.org/bio-formats/
%   - Image Processing Toolbox
%   - Deep Learning Toolbox (Level 3, matlab model only)
%   - Python + stardist/cellpose/nuclear-segmentator (Level 3, optional)
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
params.nd2_path        = 'D:\Matlab Files\Codes\SAM_Data\AnnotationData\20260413_105534_226__WellA05_Channel488,405,561,532_1_Seq0000_crop(first5FOV).nd2';

% Shared params
params.n_scan_cells    = 100;
params.n_zstack_cells  = 10;

% Channel indices (1-based: 405, 488, 561, 640 nm)
params.nucleus_channel = 1;
params.gfp_channel     = 2;
params.rfp_channel     = 3;
params.far_red_channel = 4;

% Z-stack params
params.z_planes         = 7;    % standard Z planes (L1/L2)
params.z_planes_highres = 15;   % high-res Z planes for L3 targets
params.z_step_um        = 0.5;

% Level 2 params
params.circularity_threshold = 0.85;

% ── Level 3 UNet params ───────────────────────────────────────────────────
% Model choice: 'mock' | 'stardist' | 'cellpose' | 'nuclear_seg' | 'matlab'
%
%   'mock'        No install. Rule-based proxy. Use to verify pipeline first.
%
%   'stardist'    RECOMMENDED for DAPI nuclei.
%                 pip install stardist tensorflow
%                 pyenv('Version','C:\path\to\python.exe')
%
%   'cellpose'    Good general segmentation.
%                 pip install cellpose
%                 pyenv('Version','C:\path\to\python.exe')
%
%   'nuclear_seg' Best for normality classification (direct class output).
%                 pip install nuclear-segmentator
%                 pyenv('Version','C:\path\to\python.exe')
%
%   'matlab'      Pure MATLAB. Needs Deep Learning Toolbox + trained .mat file.
%                 Set params.unet_model_path below.

params.unet_model      = 'mock';
params.unet_target     = 'round_nucleus';
params.unet_threshold  = 0.75;
params.unet_model_path = '';

% Python executable (required for stardist / cellpose / nuclear_seg)
% Uncomment and set your Python path:
% pyenv('Version', 'C:\Users\chenw\AppData\Local\Programs\Python\Python310\python.exe');

% ── Nucleus Normality Test params ─────────────────────────────────────────
% Mode: 'mock' (no model) | 'stardist' | 'cellpose' | 'nuclear_seg'
params.normality_mode            = 'mock';
params.normality_circ_threshold  = 0.80;
params.normality_solid_threshold = 0.85;
params.normality_area_ratio_hi   = 2.5;
params.normality_area_ratio_lo   = 0.30;
params.normality_model_path      = '';

% ── Ground Truth Validation params (Option 6) ─────────────────────────────
params.gt_mask_path = 'D:\Matlab Files\Codes\SAM_Data\AnnotationData\Mask\20260413_105534_226__WellA05_Channel488,405,561,532_1_Seq0000_masks_FOV1.tif';
params.gt_csv_path  = 'D:\Matlab Files\Codes\SAM_Data\AnnotationData\Meta_data\20260413_105534_226__WellA05_Channel488,405,561,532_1_Seq0000_annotations_FOV1.csv';
params.gt_nd2_path  = 'D:\Matlab Files\Codes\SAM_Data\AnnotationData\20260413_105534_226__WellA05_Channel488,405,561,532_1_Seq0000_crop(first5FOV).nd2';
params.gt_series          = 1;     % FOV1 starts at Series 1
params.gt_series_end      = 15;    % FOV1 spans Series 1-15 (15 Z planes)
params.gt_z_plane         = 1;     % not used when gt_use_mip=true
params.gt_nucleus_channel = 2;     % 405nm = plane 2 (channel order: 488,405,561,532)
params.gt_use_mip         = true;  % GT mask was generated from MIP across all Z planes

% ── 10X -> 40X Validation params (Option 7) ───────────────────────────────
params.nd2_10x_path          = 'D:\Matlab Files\Codes\SAM_Data\AnnotationData\20260412_195011_768__WellA05_Channel405,561,532_1,BF_Seq0000.nd2';
params.nd2_40x_path          = params.gt_nd2_path;
params.channel_10x           = 1;   % 405nm is 1st channel in 10X file
params.channel_40x           = 2;   % 405nm is 2nd channel in 40X file
params.n_fov_10x             = 9;
params.n_fov_40x             = 5;
params.z_planes_per_40x_fov  = 15;
script_dir           = fileparts(mfilename('fullpath'));
params.run_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
params.log_dir       = fullfile(script_dir, 'logs', params.run_timestamp);

if ~exist(params.log_dir, 'dir')
    mkdir(params.log_dir);
end
fprintf('[LOG]  Run folder: %s\n', params.log_dir);

log_file = fullfile(params.log_dir, 'run_log.txt');
diary(log_file);
fprintf('[LOG]  Console mirrored to run_log.txt\n\n');

%% ── Select Level ─────────────────────────────────────────────────────────
fprintf('\n+--------------------------------------------------+\n');
fprintf('|       SAM Benchmark Demo  --  Levels 1-3        |\n');
fprintf('+--------------------------------------------------+\n\n');
fprintf('Select benchmark level:\n');
fprintf('  [1]  Level 1 -- Random Sampling\n');
fprintf('  [2]  Level 2 -- Circularity-Guided (rule-based)\n');
fprintf('  [3]  Level 3 -- ML-Guided (UNet: %s)\n', params.unet_model);
fprintf('  [4]  Run ALL three and compare\n');
fprintf('  [5]  Nucleus Normality Test (405nm DAPI, mode: %s)\n', params.normality_mode);
fprintf('  [6]  Segmentation Validation vs Ground Truth (ROC)\n');
fprintf('  [7]  10X -> 40X GT Validation (stitched, overlap crop)\n\n');

level_choice = input('Enter choice (1/2/3/4/5/6/7): ');

%% ── Load .nd2 File (cached) ──────────────────────────────────────────────
% Option 7 handles its own loading (10X + 40X stitching)
if level_choice ~= 7
% Check if already loaded in base workspace to avoid reloading every run
if evalin('base', 'exist(''img_data_cached'',''var'')') && ...
   evalin('base', 'exist(''img_data_cached_path'',''var'')') && ...
   strcmp(evalin('base', 'img_data_cached_path'), params.nd2_path)
    fprintf('[INFO] Using cached .nd2 data (skipping reload).\n');
    fprintf('       %s\n\n', params.nd2_path);
    img_data = evalin('base', 'img_data_cached');
else
    fprintf('\n[INFO] Loading .nd2 file...\n');
    fprintf('       %s\n\n', params.nd2_path);
    img_data = load_nd2(params.nd2_path, params);
    assignin('base', 'img_data_cached',      img_data);
    assignin('base', 'img_data_cached_path', params.nd2_path);
    fprintf('[INFO] .nd2 cached for future runs.\n\n');

    % Save raw nucleus image to log (only on first load)
    nucleus_fig = figure('Visible', 'off');
    imagesc(img_data.nucleus); colormap gray; axis image off;
    title(sprintf('Nucleus channel (405nm) -- %s', params.run_timestamp), 'Interpreter', 'none');
    exportgraphics(nucleus_fig, fullfile(params.log_dir, 'nucleus_raw.png'), 'Resolution', 150);
    close(nucleus_fig);
end
end

%% ── Load GT .nd2 File for option 6 (cached separately) ──────────────────
if level_choice == 6
    if evalin('base', 'exist(''gt_img_data_cached'',''var'')') && ...
       evalin('base', 'exist(''gt_img_data_cached_path'',''var'')') && ...
       strcmp(evalin('base', 'gt_img_data_cached_path'), params.gt_nd2_path)
        fprintf('[INFO] Using cached GT .nd2 data (skipping reload).\n\n');
    else
        fprintf('[INFO] Loading GT .nd2 file...\n');
        fprintf('       %s\n\n', params.gt_nd2_path);
        gt_params                  = params;
        gt_params.nucleus_channel  = params.gt_nucleus_channel;
        gt_params.series           = params.gt_series;
        gt_params.series_end       = params.gt_series_end;
        gt_params.z_plane          = params.gt_z_plane;
        gt_params.use_mip          = params.gt_use_mip;
        gt_img_data = load_nd2(params.gt_nd2_path, gt_params);
        assignin('base', 'gt_img_data_cached',      gt_img_data);
        assignin('base', 'gt_img_data_cached_path', params.gt_nd2_path);
        fprintf('[INFO] GT .nd2 cached for future runs.\n\n');
    end
end

%% ── Run Selected Level(s) ────────────────────────────────────────────────
results_L1 = []; results_L2 = []; results_L3 = [];

switch level_choice
    case 1
        results_L1 = run_level1(img_data, params);
        display_results(results_L1, [], [], params);
        save_results(results_L1, [], [], params);

    case 2
        results_L2 = run_level2(img_data, params);
        display_results([], results_L2, [], params);
        save_results([], results_L2, [], params);

    case 3
        results_L3 = run_level3(img_data, params);
        display_results([], [], results_L3, params);
        save_results([], [], results_L3, params);

    case 4
        results_L1 = run_level1(img_data, params);
        results_L2 = run_level2(img_data, params);
        results_L3 = run_level3(img_data, params);
        display_results(results_L1, results_L2, results_L3, params);
        save_results(results_L1, results_L2, results_L3, params);

    case 5
        norm_results = run_nucleus_normality_test(img_data, params);
        save_normality_results(norm_results, params);
        fprintf('\n[DONE] Normality test complete.\n');
        fprintf('[LOG]  Outputs saved to: %s\n', params.log_dir);
        diary off;
        return;

    case 6
        % Segmentation validation vs ground truth
        gt_img_data = evalin('base', 'gt_img_data_cached');
        bench_results = run_segmentation_validation( ...
            gt_img_data, params.gt_mask_path, params.gt_csv_path, params);
        % Save
        mat_path = fullfile(params.log_dir, 'results_segmentation_validation.mat');
        save(mat_path, 'bench_results');
        fprintf('[LOG]  Saved: results_segmentation_validation.mat\n');
        fprintf('\n[DONE] Segmentation validation complete.\n');
        fprintf('[LOG]  Outputs saved to: %s\n', params.log_dir);
        diary off;
        return;

    case 7
        % 10X -> 40X GT validation using stage coordinates
        val_results = run_segmentation_validation_10x(params);
        mat_path = fullfile(params.log_dir, 'results_validation_10x.mat');
        save(mat_path, 'val_results', '-v7.3');
        fprintf('[LOG]  Saved: results_validation_10x.mat\n');
        fprintf('\n[DONE] 10X -> 40X validation complete.\n');
        fprintf('[LOG]  Outputs saved to: %s\n', params.log_dir);
        diary off;
        return;

    otherwise
        error('Invalid choice. Enter 1, 2, 3, 4, or 5.');
end

fprintf('\n[DONE] Benchmark complete.\n');
fprintf('[LOG]  Outputs saved to: %s\n', params.log_dir);
diary off;


%% ========================================================================
%  LOCAL FUNCTIONS — must all appear after the main script body
%% ========================================================================

function save_results(R1, R2, R3, params)
%SAVE_RESULTS  Save .mat and morphology CSV for each level that was run.
levels = {};
if ~isempty(R1); levels{end+1} = R1; end
if ~isempty(R2); levels{end+1} = R2; end
if ~isempty(R3); levels{end+1} = R3; end

for i = 1:numel(levels)
    R = levels{i};

    mat_path = fullfile(params.log_dir, sprintf('results_Level%d.mat', R.level));
    save(mat_path, 'R');
    fprintf('[LOG]  Saved: results_Level%d.mat\n', R.level);

    csv_path = fullfile(params.log_dir, sprintf('morphology_Level%d.csv', R.level));
    fid = fopen(csv_path, 'w');
    if R.level == 3
        fprintf(fid, 'CellID,Area_px2,Circularity,UNetScore,AspectRatio,Solidity,FocusQuality\n');
        for k = 1:numel(R.morph)
            m  = R.morph(k);
            tc = R.target_cells(k);
            fprintf(fid, '%d,%.2f,%.4f,%.4f,%.4f,%.4f,%.6f\n', ...
                    m.cell_id, m.nucleus_area, m.circularity, tc.unet_score, ...
                    m.aspect_ratio, m.solidity, m.focus_quality);
        end
    else
        fprintf(fid, 'CellID,Area_px2,Circularity,AspectRatio,Solidity,FocusQuality\n');
        for k = 1:numel(R.morph)
            m = R.morph(k);
            fprintf(fid, '%d,%.2f,%.4f,%.4f,%.4f,%.6f\n', ...
                    m.cell_id, m.nucleus_area, m.circularity, ...
                    m.aspect_ratio, m.solidity, m.focus_quality);
        end
    end
    fclose(fid);
    fprintf('[LOG]  Saved: morphology_Level%d.csv\n', R.level);
end
end


function save_normality_results(R, params)
%SAVE_NORMALITY_RESULTS  Save normality test results as .mat and CSV.
mat_path = fullfile(params.log_dir, 'results_normality.mat');
save(mat_path, 'R');
fprintf('[LOG]  Saved: results_normality.mat\n');

csv_path = fullfile(params.log_dir, 'normality_classification.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, 'CellID,Class,ClassName,Confidence,Circularity,Solidity,Area_px2,Reason\n');
for k = 1:R.n_detected
    c = R.cells(k);
    fprintf(fid, '%d,%d,%s,%.4f,%.4f,%.4f,%.0f,"%s"\n', ...
            c.id, c.norm_class, c.norm_class_name, c.norm_confidence, ...
            c.circularity, c.solidity, c.area, c.norm_reason);
end
fclose(fid);
fprintf('[LOG]  Saved: normality_classification.csv\n');
end
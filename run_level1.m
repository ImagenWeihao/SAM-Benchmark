function results = run_level1(img_data, params)
%RUN_LEVEL1  Level 1 Baseline -- Random Sampling
%
%  Workflow:
%    1. Coarse scan: segment all nuclei in the full FOV
%    2. Randomly sample N cells (params.n_scan_cells cap applied first)
%    3. From those, randomly pick M cells for Z-stack (params.n_zstack_cells)
%    4. Acquire Z-stack on selected targets
%    5. Perform cell morphology characterisation
%
%  results = run_level1(img_data, params)
fprintf('\n==============================================\n');
fprintf('  LEVEL 1 -- Random Sampling Baseline\n');
fprintf('==============================================\n');
t_start = tic;
%% -- Step 1: Coarse scan -- Otsu baseline segmentation -------------------
% L1 uses the simplest non-if-then segmentation (Otsu + connected components)
% on purpose -- it establishes a lower-tier baseline against L2 (rule-based
% circularity filter) and L3 (ML).
fprintf('\n[L1] Step 1/5 -- Coarse scan: Otsu segmentation in full FOV\n');
[all_cells, seg_mask, bw_nuclei] = segment_nuclei_l1(img_data.nucleus, params);
n_detected = numel(all_cells);
fprintf('[L1] Total nuclei detected: %d\n', n_detected);
if n_detected == 0
    error('Level 1: No nuclei detected. Check segmentation parameters.');
end
%% -- Step 2: Random subsample to n_scan_cells ----------------------------
fprintf('\n[L1] Step 2/5 -- Random subsample to %d scan candidates\n', params.n_scan_cells);
n_scan = min(params.n_scan_cells, n_detected);
scan_idx   = randperm(n_detected, n_scan);
scan_cells = all_cells(scan_idx);
fprintf('[L1] Scan pool size: %d cells\n', n_scan);
%% -- Step 3: Randomly select n_zstack_cells for Z-stack ------------------
fprintf('\n[L1] Step 3/5 -- Random selection of %d cells for Z-stack\n', params.n_zstack_cells);
n_target   = min(params.n_zstack_cells, n_scan);
target_idx = randperm(n_scan, n_target);
target_cells = scan_cells(target_idx);
fprintf('[L1] Selected cells (random): ');
fprintf('%d ', [target_cells.id]);
fprintf('\n');
%% -- Step 4: Z-stack acquisition -----------------------------------------
fprintf('\n[L1] Step 4/5 -- Z-stack acquisition on %d target cells\n', n_target);
zstack_results = acquire_zstack(target_cells, img_data, params);
%% -- Step 5: Morphology characterisation ---------------------------------
fprintf('\n[L1] Step 5/5 -- Morphology characterisation\n');
morph = characterize_morphology(target_cells, zstack_results, img_data, params);
%% -- Package results ------------------------------------------------------
elapsed = toc(t_start);
results = struct();
results.level        = 1;
results.label        = 'Level 1: Random Sampling';
results.all_cells    = all_cells;
results.seg_mask     = seg_mask;
results.bw_nuclei    = bw_nuclei;
results.scan_cells   = scan_cells;
results.target_cells = target_cells;
results.zstack       = zstack_results;
results.morph        = morph;
results.n_detected   = n_detected;
results.n_scanned    = n_scan;
results.n_targeted   = n_target;
results.elapsed_sec  = elapsed;
% Selection stats
circ_all      = [all_cells.circularity];
circ_selected = [target_cells.circularity];
results.stats.circ_mean_all      = mean(circ_all);
results.stats.circ_mean_selected = mean(circ_selected);
results.stats.circ_std_selected  = std(circ_selected);
fprintf('\n[L1] -- Summary -----------------------------\n');
fprintf('[L1]   Detected:  %d nuclei\n',   n_detected);
fprintf('[L1]   Scanned:   %d (random)\n', n_scan);
fprintf('[L1]   Targeted:  %d (random)\n', n_target);
fprintf('[L1]   Mean circularity (all):      %.3f\n', results.stats.circ_mean_all);
fprintf('[L1]   Mean circularity (selected): %.3f +/- %.3f\n', ...
        results.stats.circ_mean_selected, results.stats.circ_std_selected);
fprintf('[L1]   Elapsed: %.2f s\n', elapsed);
fprintf('[L1] ----------------------------------------\n');
end
%% ======================================================================
function [cells, seg_mask, bw_nuclei] = segment_nuclei_l1(nucleus_img, params) %#ok<INUSD>
%SEGMENT_NUCLEI_L1  Otsu + connected-components baseline (Level 1).
%  Intentionally minimal: no morphological cleanup, no area filter, no
%  watershed splitting, no circularity rule.
fprintf('[SEG-L1] Otsu + connected components...\n');
level     = graythresh(nucleus_img);
bw_nuclei = imbinarize(nucleus_img, level);
seg_mask  = uint16(bwlabel(bw_nuclei));
props     = regionprops(seg_mask, nucleus_img, ...
    'Area', 'Centroid', 'Perimeter', 'BoundingBox', 'MeanIntensity', 'Solidity');
n_detected = numel(props);
fprintf('[SEG-L1] Detected %d connected components.\n', n_detected);
if n_detected == 0
    cells = struct([]);
    return;
end
cells = struct( ...
    'id',          num2cell((1:n_detected)'), ...
    'centroid',    {props.Centroid}', ...
    'area',        num2cell([props.Area]'), ...
    'perimeter',   num2cell([props.Perimeter]'), ...
    'circularity', num2cell(zeros(n_detected,1)), ...
    'solidity',    num2cell([props.Solidity]'), ...
    'bbox',        {props.BoundingBox}', ...
    'mean_int',    num2cell([props.MeanIntensity]'));
for i = 1:n_detected
    p = cells(i).perimeter;
    a = cells(i).area;
    if p > 0
        cells(i).circularity = min(1, (4 * pi * a) / (p^2));
    else
        cells(i).circularity = 0;
    end
end
end

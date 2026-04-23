function results = run_level3(img_data, params)
%RUN_LEVEL3  Level 3 -- ML-Guided Target Selection via UNet Classifier
%
%  Workflow:
%    1. Coarse scan: segment all nuclei in full FOV (same as L2; L1 uses Otsu)
%    2. Run UNet classifier on each candidate cell (user-selected model)
%    3. Rank candidates by classification confidence score
%    4. Select top-N cells for high-resolution Z-stack
%    5. Morphology characterisation (shared)
%
%  The UNet model and target morphology class are set in params:
%    params.unet_model      -- 'cellpose' | 'stardist' | 'matlab' | 'mock'
%    params.unet_target     -- string describing target class, e.g. 'mitotic'
%    params.unet_threshold  -- minimum confidence score to qualify (0-1)
%
%  results = run_level3(img_data, params)
fprintf('\n==============================================\n');
fprintf('  LEVEL 3 -- ML-Guided Target Selection\n');
fprintf('  Model:  %s\n', params.unet_model);
fprintf('  Target: %s\n', params.unet_target);
fprintf('  Confidence threshold: %.2f\n', params.unet_threshold);
fprintf('==============================================\n');
t_start = tic;
%% -- Step 1: Coarse scan -- detect all nuclei -----------------------------
fprintf('\n[L3] Step 1/5 -- Coarse scan: detecting nuclei in full FOV\n');
[all_cells, seg_mask, bw_nuclei] = segment_nuclei(img_data.nucleus, params);
n_detected = numel(all_cells);
fprintf('[L3] Total nuclei detected: %d\n', n_detected);
if n_detected == 0
    error('Level 3: No nuclei detected. Check segmentation parameters.');
end
%% -- Step 2: Subsample to scan pool (same budget as L1/L2) ---------------
fprintf('\n[L3] Step 2/5 -- Subsample to %d scan candidates\n', params.n_scan_cells);
n_scan     = min(params.n_scan_cells, n_detected);
scan_idx   = randperm(n_detected, n_scan);
scan_cells = all_cells(scan_idx);
fprintf('[L3] Scan pool: %d cells\n', n_scan);
%% -- Step 3: UNet classification -----------------------------------------
fprintf('\n[L3] Step 3/5 -- UNet classification (%s) on scan pool\n', params.unet_model);
[scan_cells, model_info] = classify_cells_unet(scan_cells, img_data, params);
% Report score distribution
scores = [scan_cells.unet_score];
fprintf('[L3] Confidence scores -- min: %.3f  max: %.3f  mean: %.3f\n', ...
        min(scores), max(scores), mean(scores));
%% -- Step 4: Select top-N by confidence score ----------------------------
fprintf('\n[L3] Step 4/5 -- Selecting targets (confidence >= %.2f)\n', params.unet_threshold);
pass_mask  = scores >= params.unet_threshold;
pass_cells = scan_cells(pass_mask);
n_pass     = sum(pass_mask);
fprintf('[L3] Cells passing threshold: %d / %d\n', n_pass, n_scan);
if n_pass > 0
    [~, sort_ord] = sort([pass_cells.unet_score], 'descend');
    pass_cells     = pass_cells(sort_ord);
    n_target       = min(params.n_zstack_cells, n_pass);
    target_cells   = pass_cells(1:n_target);
    selection_mode = 'threshold + top-N by confidence';
else
    warning('[L3] No cells exceed confidence threshold %.2f. Using top-%d by score.', ...
            params.unet_threshold, params.n_zstack_cells);
    [~, sort_ord] = sort(scores, 'descend');
    n_target     = min(params.n_zstack_cells, n_scan);
    target_cells = scan_cells(sort_ord(1:n_target));
    selection_mode = 'fallback top-N by confidence';
end
fprintf('[L3] Selected %d cells | mode: %s\n', n_target, selection_mode);
fprintf('[L3] Confidence values: ');
fprintf('%.3f ', [target_cells.unet_score]);
fprintf('\n');
%% -- Step 5: Z-stack + morphology (shared) -------------------------------
fprintf('\n[L3] Step 5a/5 -- High-res Z-stack acquisition\n');
% Level 3 uses more Z planes for high-res imaging
params_hr          = params;
params_hr.z_planes = params.z_planes_highres;
zstack_results     = acquire_zstack(target_cells, img_data, params_hr);
fprintf('\n[L3] Step 5b/5 -- Morphology characterisation\n');
morph = characterize_morphology(target_cells, zstack_results, img_data, params);
%% -- Package results ------------------------------------------------------
elapsed = toc(t_start);
results                    = struct();
results.level              = 3;
results.label              = sprintf('Level 3: ML-Guided (%s)', params.unet_model);
results.all_cells          = all_cells;
results.seg_mask           = seg_mask;
results.bw_nuclei          = bw_nuclei;
results.scan_cells         = scan_cells;
results.target_cells       = target_cells;
results.zstack             = zstack_results;
results.morph              = morph;
results.model_info         = model_info;
results.n_detected         = n_detected;
results.n_scanned          = n_scan;
results.n_passed_threshold = n_pass;
results.n_targeted         = n_target;
results.selection_mode     = selection_mode;
results.elapsed_sec        = elapsed;
% Stats
scores_all      = [all_cells.circularity];
scores_selected = [target_cells.unet_score];
circ_selected   = [target_cells.circularity];
results.stats.unet_score_mean    = mean(scores_selected);
results.stats.unet_score_std     = std(scores_selected);
results.stats.circ_mean_all      = mean(scores_all);
results.stats.circ_mean_selected = mean(circ_selected);
results.stats.pct_above_threshold= 100 * n_pass / n_scan;
fprintf('\n[L3] -- Summary ----------------------------------\n');
fprintf('[L3]   Detected:           %d nuclei\n',   n_detected);
fprintf('[L3]   Scanned:            %d\n',           n_scan);
fprintf('[L3]   Passed threshold:   %d (%.1f%%)\n', n_pass, results.stats.pct_above_threshold);
fprintf('[L3]   Targeted:           %d\n',           n_target);
fprintf('[L3]   Mean UNet score (selected): %.3f +/- %.3f\n', ...
        results.stats.unet_score_mean, results.stats.unet_score_std);
fprintf('[L3]   Mean circularity (selected): %.3f\n', results.stats.circ_mean_selected);
fprintf('[L3]   Model: %s  |  Target: %s\n', params.unet_model, params.unet_target);
fprintf('[L3]   Elapsed: %.2f s\n', elapsed);
fprintf('[L3] -------------------------------------------------\n');
end

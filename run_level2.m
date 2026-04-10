function results = run_level2(img_data, params)
%RUN_LEVEL2  Level 2 — Circularity-Guided Selection (If-Then Rule)
%
%  Workflow:
%    1. Coarse scan: segment all nuclei in the full FOV
%    2. Subsample to n_scan_cells (same pool as Level 1 for fair comparison)
%    3. Apply circularity threshold rule to identify target cells
%       — Selects cells with circularity >= params.circularity_threshold
%       — If fewer than n_zstack_cells pass, top-N by circularity are used
%       — If more pass than n_zstack_cells, top-N by circularity are taken
%    4. Acquire Z-stack on selected targets
%    5. Perform cell morphology characterisation
%
%  results = run_level2(img_data, params)

fprintf('\n══════════════════════════════════════════════\n');
fprintf('  LEVEL 2 — Circularity-Guided Selection\n');
fprintf('  Threshold: circularity >= %.2f\n', params.circularity_threshold);
fprintf('══════════════════════════════════════════════\n');
t_start = tic;

%% ── Step 1: Coarse scan — detect all nuclei ─────────────────────────────
fprintf('\n[L2] Step 1/5 — Coarse scan: detecting nuclei in full FOV\n');
[all_cells, seg_mask, bw_nuclei] = segment_nuclei(img_data.nucleus, params);

n_detected = numel(all_cells);
fprintf('[L2] Total nuclei detected: %d\n', n_detected);

if n_detected == 0
    error('Level 2: No nuclei detected. Check segmentation parameters.');
end

%% ── Step 2: Subsample to n_scan_cells (same budget as Level 1) ──────────
fprintf('\n[L2] Step 2/5 — Random subsample to %d scan candidates\n', params.n_scan_cells);

n_scan = min(params.n_scan_cells, n_detected);
scan_idx   = randperm(n_detected, n_scan);
scan_cells = all_cells(scan_idx);

fprintf('[L2] Scan pool: %d cells\n', n_scan);

%% ── Step 3: If-then circularity selection ───────────────────────────────
fprintf('\n[L2] Step 3/5 — Applying circularity rule (threshold = %.2f)\n', ...
        params.circularity_threshold);

circ_vals   = [scan_cells.circularity];
pass_mask   = circ_vals >= params.circularity_threshold;
pass_cells  = scan_cells(pass_mask);
n_pass      = sum(pass_mask);

fprintf('[L2] Cells passing circularity threshold: %d / %d\n', n_pass, n_scan);

% Sort passing cells by circularity (descending) and take top-N
if n_pass > 0
    [~, sort_ord] = sort([pass_cells.circularity], 'descend');
    pass_cells = pass_cells(sort_ord);
    n_target   = min(params.n_zstack_cells, n_pass);
    target_cells = pass_cells(1:n_target);
    selection_mode = 'threshold + top-N';
else
    % Fallback: no cells pass threshold — use top-N by circularity from scan pool
    warning('[L2] No cells exceed circularity threshold %.2f. Falling back to top-%d by circularity.', ...
            params.circularity_threshold, params.n_zstack_cells);
    [~, sort_ord] = sort(circ_vals, 'descend');
    n_target   = min(params.n_zstack_cells, n_scan);
    target_cells = scan_cells(sort_ord(1:n_target));
    selection_mode = 'fallback top-N (no threshold passed)';
end

fprintf('[L2] Selection mode: %s\n', selection_mode);
fprintf('[L2] Selected cells (IDs): ');
fprintf('%d ', [target_cells.id]);
fprintf('\n');
fprintf('[L2] Selected circularity values: ');
fprintf('%.3f ', [target_cells.circularity]);
fprintf('\n');

%% ── Step 4: Z-stack acquisition ─────────────────────────────────────────
fprintf('\n[L2] Step 4/5 — Z-stack acquisition on %d target cells\n', n_target);
zstack_results = acquire_zstack(target_cells, img_data, params);

%% ── Step 5: Morphology characterisation ─────────────────────────────────
fprintf('\n[L2] Step 5/5 — Morphology characterisation\n');
morph = characterize_morphology(target_cells, zstack_results, img_data, params);

%% ── Package results ──────────────────────────────────────────────────────
elapsed = toc(t_start);

results = struct();
results.level          = 2;
results.label          = 'Level 2: Circularity-Guided';
results.all_cells      = all_cells;
results.seg_mask       = seg_mask;
results.bw_nuclei      = bw_nuclei;
results.scan_cells     = scan_cells;
results.target_cells   = target_cells;
results.zstack         = zstack_results;
results.morph          = morph;
results.n_detected     = n_detected;
results.n_scanned      = n_scan;
results.n_targeted     = n_target;
results.n_passed_threshold = n_pass;
results.selection_mode = selection_mode;
results.elapsed_sec    = elapsed;

circ_all      = [all_cells.circularity];
circ_selected = [target_cells.circularity];

results.stats.circ_mean_all      = mean(circ_all);
results.stats.circ_mean_selected = mean(circ_selected);
results.stats.circ_std_selected  = std(circ_selected);
results.stats.circ_threshold     = params.circularity_threshold;
results.stats.pct_above_threshold= 100 * n_pass / n_scan;

fprintf('\n[L2] ── Summary ─────────────────────────────\n');
fprintf('[L2]   Detected:  %d nuclei\n',      n_detected);
fprintf('[L2]   Scanned:   %d (subsample)\n', n_scan);
fprintf('[L2]   Passed threshold (≥%.2f): %d (%.1f%%)\n', ...
        params.circularity_threshold, n_pass, results.stats.pct_above_threshold);
fprintf('[L2]   Targeted:  %d\n',                  n_target);
fprintf('[L2]   Mean circularity (all):      %.3f\n', results.stats.circ_mean_all);
fprintf('[L2]   Mean circularity (selected): %.3f ± %.3f\n', ...
        results.stats.circ_mean_selected, results.stats.circ_std_selected);
fprintf('[L2]   Elapsed: %.2f s\n', elapsed);
fprintf('[L2] ────────────────────────────────────────\n');
end

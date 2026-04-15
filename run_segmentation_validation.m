function results = run_segmentation_validation(img_data, mask_path, csv_path, params)
%RUN_SEGMENTATION_VALIDATION  Compare segmentation methods vs annotated ground truth.
%
%  img_data   — pre-loaded image struct from load_nd2 (cached in run_benchmark)
%  mask_path  — multi-frame TIFF: Frame 0 = nucleus label mask (uint16)
%  csv_path   — annotation CSV

fprintf('\n==============================================\n');
fprintf('  SEGMENTATION BENCHMARK vs GROUND TRUTH\n');
fprintf('  FOV1  |  405nm DAPI channel  |  Series=%d\n', params.gt_series);
fprintf('==============================================\n\n');

%% ── Step 1: Load ground truth mask and annotations ───────────────────────
fprintf('[BENCH] Loading ground truth...\n');

% Load nucleus mask (Frame 0)
gt_mask = load_gt_mask(mask_path, 0);
fprintf('[BENCH] GT mask: %dx%d, %d nuclei labelled\n', ...
        size(gt_mask,1), size(gt_mask,2), max(gt_mask(:)));

% Load annotations CSV
annot = readtable(csv_path);
fprintf('[BENCH] Annotations: %d cells\n', height(annot));
fprintf('[BENCH] Morphology breakdown:\n');
morphs = unique(annot.morphology);
for m = 1:numel(morphs)
    n_m = sum(strcmp(annot.morphology, morphs{m}));
    fprintf('          %-12s : %d\n', morphs{m}, n_m);
end

%% ── Step 2: Use pre-loaded image data ───────────────────────────────────
fprintf('[BENCH] Using pre-loaded image (W=%d H=%d C=%d)\n', ...
        img_data.metadata.img_width, img_data.metadata.img_height, img_data.n_channels);
nucleus_img = img_data.nucleus;

% Resize nucleus image to match GT mask if needed
[Hm, Wm] = size(gt_mask);
[Hi, Wi] = size(nucleus_img);
if Hi ~= Hm || Wi ~= Wm
    fprintf('[BENCH] Resizing image from %dx%d to %dx%d to match GT mask\n', ...
            Wi, Hi, Wm, Hm);
    nucleus_img = imresize(nucleus_img, [Hm Wm]);
end

%% ── Step 3: Run each segmentation method ────────────────────────────────
fprintf('\n[BENCH] Running segmentation methods...\n');

% Method 1: MATLAB adaptive threshold (used by L1/L2)
fprintf('[BENCH] Method 1: MATLAB segmentation (L1/L2 baseline)...\n');
[cells_matlab, seg_matlab] = segment_nuclei(nucleus_img, params);
fprintf('[BENCH]   Detected %d nuclei\n', numel(cells_matlab));

% Method 2: Circularity-filtered (L2)
% Use lower threshold for this dataset (mean circ ~0.635)
circ_thresh_val = min(params.circularity_threshold, 0.60);
fprintf('[BENCH] Method 2: Circularity-filtered (L2, thresh=%.2f)...\n', circ_thresh_val);
circ_vals    = [cells_matlab.circularity];
pass_mask_L2 = circ_vals >= circ_thresh_val;
seg_L2       = zeros(size(seg_matlab), 'uint16');
cell_ids_L2  = find(pass_mask_L2);
for k = 1:numel(cell_ids_L2)
    seg_L2(seg_matlab == cell_ids_L2(k)) = k;
end
fprintf('[BENCH]   %d/%d cells pass circularity threshold\n', ...
        sum(pass_mask_L2), numel(cells_matlab));

% Method 3: StarDist via bridge
fprintf('[BENCH] Method 3: StarDist (L3)...\n');
seg_stardist = run_stardist_full_image(nucleus_img, params);
n_sd = max(seg_stardist(:));
fprintf('[BENCH]   StarDist detected %d nuclei\n', n_sd);

%% ── Step 4: Compute IoU-based ROC for each method ───────────────────────
fprintf('\n[BENCH] Computing ROC curves...\n');

methods = {'MATLAB (L1/L2)', 'Circularity-filtered (L2)', 'StarDist (L3)'};
segs    = {seg_matlab, seg_L2, seg_stardist};
colors  = {[0.3 0.6 1.0], [1.0 0.4 0.3], [0.2 0.8 0.3]};
roc     = cell(1, 3);

for m = 1:3
    roc{m} = compute_roc(gt_mask, segs{m}, annot);
    fprintf('[BENCH] %s — AUC = %.3f  |  F1@0.5 = %.3f\n', ...
            methods{m}, roc{m}.auc, roc{m}.f1_at_half);
end

%% ── Step 5: Display results ──────────────────────────────────────────────
plot_benchmark_results(nucleus_img, gt_mask, segs, roc, methods, colors, annot, params);

%% ── Package output ───────────────────────────────────────────────────────
results.roc        = roc;
results.methods    = methods;
results.gt_mask    = gt_mask;
results.segs       = segs;
results.annot      = annot;
results.nucleus_img= nucleus_img;
end


%% ======================================================================
function gt_mask = load_gt_mask(mask_path, frame_idx)
%LOAD_GT_MASK  Load a specific frame from a multi-frame TIFF.
% Uses imread with frame index. Frame 0 = nucleus mask.
info     = imfinfo(mask_path);
gt_mask  = imread(mask_path, frame_idx + 1);   % imread is 1-based
gt_mask  = uint16(gt_mask);
end


%% ======================================================================
function seg_out = run_stardist_full_image(nucleus_img, params)
%RUN_STARDIST_FULL_IMAGE  Run StarDist using flat array transfer to fix MATLAB/Python memory order.

seg_out = zeros(size(nucleus_img), 'uint16');

bridge_path = fullfile(fileparts(mfilename('fullpath')), 'stardist_bridge.py');
if ~isfile(bridge_path)
    warning('[StarDist] Bridge not found at: %s', bridge_path);
    return;
end

try
    spec   = py.importlib.util.spec_from_file_location('stardist_bridge', bridge_path);
    bridge = py.importlib.util.module_from_spec(spec);
    spec.loader.exec_module(bridge);
    load_fn    = py.getattr(bridge, 'load_model');
    model      = feval(load_fn, '2D_versatile_fluo');
    predict_fn = py.getattr(bridge, 'predict_full');
    np         = py.importlib.import_module('numpy');

    % Pass as flat column-major vector + explicit dims to preserve MATLAB memory order
    img_f32  = single(nucleus_img) ./ (max(single(nucleus_img(:))) + eps);
    [H, W]   = size(img_f32);
    py_flat  = np.array(img_f32(:));   % MATLAB flattens column-major by default

    result   = feval(predict_fn, model, py_flat, int32(H), int32(W));

    % Read flat labels and reshape back to 2D
    labels_list = result{'labels_flat'};
    labels_flat = double(py.array.array('l', labels_list));
    seg_out     = uint16(reshape(labels_flat, H, W));

catch ME
    warning('[StarDist] Full image prediction failed: %s', ME.message);
end
end


%% ======================================================================
function roc = compute_roc(gt_mask, pred_mask, annot)
%COMPUTE_ROC  Compute ROC curve using per-nucleus IoU matching.
%
%  For each IoU threshold t in [0,1]:
%    TP = GT nuclei matched by a prediction with IoU >= t
%    FN = GT nuclei not matched
%    FP = predicted regions not matching any GT nucleus
%    TN = background pixels correctly not predicted
%  TPR = TP / (TP + FN)
%  FPR = FP / (FP + TN)

thresholds = 0:0.02:1.0;
n_thresh   = numel(thresholds);
tpr        = zeros(1, n_thresh);
fpr        = zeros(1, n_thresh);
precision  = zeros(1, n_thresh);

gt_ids   = unique(gt_mask(gt_mask > 0));
pred_ids = unique(pred_mask(pred_mask > 0));
n_gt     = numel(gt_ids);
n_pred   = numel(pred_ids);

% Pre-compute IoU matrix [n_gt x n_pred]
fprintf('    Computing IoU matrix (%d GT x %d pred)...\n', n_gt, n_pred);
iou_matrix = zeros(n_gt, n_pred);
for gi = 1:n_gt
    gt_mask_i = (gt_mask == gt_ids(gi));
    for pi = 1:n_pred
        pred_mask_p = (pred_mask == pred_ids(pi));
        intersection = sum(gt_mask_i(:) & pred_mask_p(:));
        if intersection > 0
            union = sum(gt_mask_i(:) | pred_mask_p(:));
            iou_matrix(gi, pi) = intersection / union;
        end
    end
end

% Total background pixels for FPR calculation
total_bg = sum(gt_mask(:) == 0);

for ti = 1:n_thresh
    t  = thresholds(ti);
    matched_gt   = false(n_gt, 1);
    matched_pred = false(n_pred, 1);

    % Greedy matching: highest IoU first
    [sorted_iou, sort_idx] = sort(iou_matrix(:), 'descend');
    for si = 1:numel(sorted_iou)
        if sorted_iou(si) < t; break; end
        [gi, pi] = ind2sub([n_gt, n_pred], sort_idx(si));
        if ~matched_gt(gi) && ~matched_pred(pi)
            matched_gt(gi)   = true;
            matched_pred(pi) = true;
        end
    end

    TP = sum(matched_gt);
    FN = n_gt - TP;
    FP = sum(~matched_pred);

    % FPR: false positive pixels / total background pixels
    fp_px = 0;
    for pi = 1:n_pred
        if ~matched_pred(pi)
            fp_px = fp_px + sum(pred_mask(:) == pred_ids(pi));
        end
    end

    tpr(ti)       = TP / (n_gt + eps);
    fpr(ti)       = fp_px / (total_bg + eps);
    if TP + FP > 0
        precision(ti) = TP / (TP + FP);
    end
end

% AUC via trapezoidal integration
[fpr_s, idx] = sort(fpr);
tpr_s        = tpr(idx);
auc          = trapz(fpr_s, tpr_s);

% F1 at IoU=0.5
t50_idx      = find(thresholds >= 0.5, 1);
f1_at_half   = 2 * precision(t50_idx) * tpr(t50_idx) / ...
               (precision(t50_idx) + tpr(t50_idx) + eps);

roc.thresholds = thresholds;
roc.tpr        = tpr;
roc.fpr        = fpr;
roc.precision  = precision;
roc.auc        = auc;
roc.f1_at_half = f1_at_half;
roc.n_gt       = n_gt;
roc.iou_matrix = iou_matrix;
end


%% ======================================================================
function plot_benchmark_results(nucleus_img, gt_mask, segs, roc, methods, colors, annot, params)
%PLOT_BENCHMARK_RESULTS  5-panel figure: overlays + ROC curve.

fig = figure('Name', 'Segmentation Benchmark vs Ground Truth', ...
             'NumberTitle', 'off', 'Units', 'normalized', ...
             'Position', [0.02 0.02 0.96 0.90]);
sgtitle('Segmentation Validation — FOV1  |  405nm DAPI  |  Methods vs Ground Truth', ...
        'FontSize', 14, 'FontWeight', 'bold');

method_short = {'MATLAB (L1/L2)', 'Circ.-filtered (L2)', 'StarDist (L3)'};

%% Panels 1-3: Segmentation overlays vs GT
for m = 1:3
    ax = subplot(2, 4, m);

    % Draw predicted boundaries in method color over nucleus image
    pred_bw   = segs{m} > 0;
    pred_edge = edge(pred_bw, 'canny');
    gt_bw     = gt_mask > 0;
    gt_edge   = edge(gt_bw, 'canny');

    img_data_rgb = repmat(mat2gray(nucleus_img), [1 1 3]);
    img_data_rgb(:,:,1) = img_data_rgb(:,:,1) + pred_edge * colors{m}(1) * 0.8;
    img_data_rgb(:,:,2) = img_data_rgb(:,:,2) + pred_edge * colors{m}(2) * 0.8;
    img_data_rgb(:,:,3) = img_data_rgb(:,:,3) + pred_edge * colors{m}(3) * 0.8;
    img_data_rgb(:,:,:) = img_data_rgb(:,:,:) + repmat(gt_edge * 0.5, [1 1 3]);
    imagesc(ax, min(1, img_data_rgb)); axis(ax,'image'); axis(ax,'off');
    hold(ax, 'on');

    % Determine which GT nuclei were matched at IoU=0.5
    iou_mat  = roc{m}.iou_matrix;
    gt_ids   = unique(gt_mask(gt_mask > 0));
    n_gt     = numel(gt_ids);
    matched  = false(n_gt, 1);
    pred_ids_all = unique(segs{m}(segs{m} > 0));
    n_pred   = numel(pred_ids_all);
    if n_pred > 0 && size(iou_mat,1) == n_gt && size(iou_mat,2) == n_pred
        [sorted_iou, sort_idx] = sort(iou_mat(:), 'descend');
        matched_pred = false(n_pred, 1);
        for si = 1:numel(sorted_iou)
            if sorted_iou(si) < 0.5; break; end
            [gi, pi] = ind2sub([n_gt, n_pred], sort_idx(si));
            if ~matched(gi) && ~matched_pred(pi)
                matched(gi) = true;
                matched_pred(pi) = true;
            end
        end
    end

    % Draw all GT centroids — green tick for detected, red X for missed
    for k = 1:height(annot)
        cx = annot.centroid_x(k);
        cy = annot.centroid_y(k);
        % Match annotation cell_id to GT mask index
        gt_row = find(gt_ids == annot.cell_id(k), 1);
        if ~isempty(gt_row) && matched(gt_row)
            % Detected — small white cross
            plot(ax, cx, cy, 'w+', 'MarkerSize', 6, 'LineWidth', 1);
        else
            % Missed — red X with cell ID label
            plot(ax, cx, cy, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
            text(ax, cx+12, cy, num2str(annot.cell_id(k)), ...
                 'Color', [1 0.3 0.3], 'FontSize', 7, 'FontWeight', 'bold');
        end
    end

    n_missed = sum(~matched(1:min(n_gt, height(annot))));
    title(ax, {method_short{m}, sprintf('AUC=%.3f  F1@0.5=%.3f  missed=%d', ...
          roc{m}.auc, roc{m}.f1_at_half, n_missed)}, 'FontSize', 9);
end

%% Panel 4: GT mask with morphology labels
ax4 = subplot(2, 4, 4);
imagesc(ax4, nucleus_img); colormap(ax4, 'gray'); axis(ax4,'image'); axis(ax4,'off');
hold(ax4, 'on');
morph_colors = containers.Map(...
    {'normal','blebbing','micronuclei','unsure'}, ...
    {[0.2 0.9 0.2], [1.0 0.4 0.0], [0.9 0.1 0.9], [0.7 0.7 0.7]});
for k = 1:height(annot)
    cx  = annot.centroid_x(k);
    cy  = annot.centroid_y(k);
    mor = annot.morphology{k};
    if isKey(morph_colors, mor)
        col = morph_colors(mor);
    else
        col = [0.7 0.7 0.7];
    end
    plot(ax4, cx, cy, 'o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 8);
    text(ax4, cx+10, cy, num2str(annot.cell_id(k)), 'Color','white', 'FontSize', 6);
end
title(ax4, {'Ground Truth Annotations', 'green=normal  orange=bleb  pink=micro  gray=unsure'}, ...
      'FontSize', 8);

%% Panel 5 (bottom, spanning cols 1-2): ROC curves
ax5 = subplot(2, 4, [5 6]);
hold(ax5, 'on');
% Diagonal reference line
plot(ax5, [0 1], [0 1], 'k--', 'LineWidth', 1, 'DisplayName', 'Random');
for m = 1:3
    [fpr_s, idx] = sort(roc{m}.fpr);
    tpr_s        = roc{m}.tpr(idx);
    plot(ax5, fpr_s, tpr_s, '-', 'Color', colors{m}, 'LineWidth', 2.5, ...
         'DisplayName', sprintf('%s (AUC=%.3f)', method_short{m}, roc{m}.auc));
    % Mark F1@0.5 point
    t50 = find(roc{m}.thresholds >= 0.5, 1);
    plot(ax5, roc{m}.fpr(t50), roc{m}.tpr(t50), 'o', ...
         'Color', colors{m}, 'MarkerSize', 10, 'MarkerFaceColor', colors{m}, ...
         'HandleVisibility', 'off');
end
xlabel(ax5, 'False Positive Rate'); ylabel(ax5, 'True Positive Rate');
legend(ax5, 'show', 'Location', 'southeast', 'FontSize', 9);
title(ax5, {'ROC Curves — Nucleus Detection', '(dots = operating point at IoU threshold 0.5)'}, ...
      'FontSize', 10);
xlim(ax5, [0 1]); ylim(ax5, [0 1]);
grid(ax5, 'on'); axis(ax5, 'square');

%% Panel 6 (bottom, spanning cols 3-4): Per-morphology detection bar chart
ax6 = subplot(2, 4, [7 8]);
morph_list = {'normal', 'blebbing', 'micronuclei', 'unsure'};
bar_data   = zeros(3, numel(morph_list));

for m = 1:3
    for mi = 1:numel(morph_list)
        mor_idx = strcmp(annot.morphology, morph_list{mi});
        if sum(mor_idx) == 0; continue; end
        % Check how many of this morphology type were detected at IoU>=0.5
        gt_ids_mor = annot.cell_id(mor_idx);
        detected   = 0;
        for gi = 1:numel(gt_ids_mor)
            gt_row = find(unique(roc{m}.iou_matrix ~= 0));
            if gi <= size(roc{m}.iou_matrix, 1)
                if max(roc{m}.iou_matrix(gi,:)) >= 0.5
                    detected = detected + 1;
                end
            end
        end
        bar_data(m, mi) = detected / sum(mor_idx) * 100;
    end
end

b = bar(ax6, bar_data', 'grouped');
for m = 1:3; b(m).FaceColor = colors{m}; end
set(ax6, 'XTickLabel', morph_list, 'XTick', 1:numel(morph_list));
ylabel(ax6, 'Detection Rate (%)');
legend(ax6, method_short, 'Location', 'best', 'FontSize', 8);
title(ax6, 'Detection Rate by Morphology Class (IoU >= 0.5)', 'FontSize', 10);
ylim(ax6, [0 110]); grid(ax6, 'on');

drawnow;

% Save
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_path = fullfile(params.log_dir, ...
        sprintf('segmentation_benchmark_%s.png', params.run_timestamp));
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Benchmark figure saved.\n');
end
end
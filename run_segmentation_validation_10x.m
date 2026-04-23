function results = run_segmentation_validation_10x(params)
%RUN_SEGMENTATION_VALIDATION_10X  Validate 10X segmentation using 40X GT mask.
%
%  Workflow:
%    1. Load 10X .nd2 and stitch all FOVs into a single mosaic using stage
%       coordinates + cameraTransformationMatrix from metadata.
%    2. Load 40X .nd2 FOVs, compute MIP per FOV, and map each 40X FOV to its
%       physical location in the 10X mosaic canvas.
%    3. Place 40X GT mask into the 10X canvas (downsampled to 10X pixel scale).
%    4. Crop the overlap region: area where BOTH 10X image AND 40X GT exist.
%    5. Run L1 (MATLAB seg), L2 (circularity), L3 (StarDist) on the 10X overlap.
%    6. Compute IoU-based ROC against the 40X GT mask (same overlap region).
%    7. Display diagnostic + results.
%
%  params fields used:
%    .nd2_10x_path            path to 10X .nd2
%    .nd2_40x_path            path to 40X .nd2
%    .gt_mask_path            path to 40X GT mask tif (Frame 0 = nuclei)
%    .gt_csv_path             path to annotation CSV (optional)
%    .channel_10x             1-based DAPI channel index in 10X file
%    .channel_40x             1-based DAPI channel index in 40X file
%    .n_fov_10x               number of 10X FOVs (9)
%    .n_fov_40x               number of 40X FOVs (5, each with 15 Z planes)
%    .z_planes_per_40x_fov    15
%    .gt_roi_enabled          [true] draw interactive ROI before accuracy scoring
%                             after NCC refinement, before segmentation/scoring.
%                             Draw a rectangle on the 10X image, double-click to
%                             confirm.  Only the selected sub-region is scored,
%                             eliminating false-misses from 10X FOV edges that
%                             fall outside the 40X acquisition coverage.
if nargin < 1; error('params struct required'); end
fprintf('\n====================================================\n');
fprintf('  SEGMENTATION VALIDATION  10X input vs 40X GT mask\n');
fprintf('====================================================\n\n');
%% -- Step 1: Load 10X and build stitched mosaic --------------------------
fprintf('[10X]  Loading overview and building mosaic...\n');
[mosaic_10x, canvas_info, px_10x] = build_10x_mosaic(params);
fprintf('[10X]  Mosaic: %d x %d px  (%.0f x %.0f um)\n', ...
        size(mosaic_10x,2), size(mosaic_10x,1), ...
        size(mosaic_10x,2)*px_10x, size(mosaic_10x,1)*px_10x);
%% -- Step 2: Load 40X FOVs, get stage coords, compute MIPs ---------------
fprintf('\n[40X]  Loading FOVs + computing MIPs + reading stage coords...\n');
[fov_40x_mips, fov_40x_coords, px_40x] = load_40x_mips(params);
%% -- Step 3: Place 40X GT mask into 10X canvas ---------------------------
% The 40X FOV can extend beyond the stitched 10X canvas (the 40X was
% imaged at a physical location whose full 333x332 um footprint partially
% falls outside the 10X tile coverage).  place_gt_in_canvas now pads the
% mosaic and canvas as needed so the full 40X mask -- every labelled cell
% -- is placed without clipping.  Padded regions are black in the 10X
% image but keep their GT labels, so the user sees the entire 40X FOV.
fprintf('\n[GT]   Placing 40X GT mask into 10X canvas...\n');
[gt_mask_canvas, fov_rect_canvas, mosaic_10x, canvas_info] = ...
    place_gt_in_canvas( ...
        params, fov_40x_coords, canvas_info, mosaic_10x, px_10x, px_40x);
%% -- Step 4: Crop overlap region -----------------------------------------
% Crop the FULL 40X FOV rectangle (not just the bbox of GT-labelled cells)
% intersected with the 10X mosaic's imaged region.  This reveals the entire
% physical area covered by the 40X objective, including any 10X nuclei that
% fall inside the 40X FOV but were not annotated in the GT mask -- the user
% can see them and judge coverage/alignment directly.
fprintf('\n[CROP] Cropping overlap region (full 40X FOV rectangle)...\n');
[crop_10x, crop_gt, crop_bbox] = crop_overlap( ...
    mosaic_10x, gt_mask_canvas, fov_rect_canvas);
fprintf('[CROP] Overlap region: %d x %d px (%.0f x %.0f um)\n', ...
        size(crop_10x,2), size(crop_10x,1), ...
        size(crop_10x,2)*px_10x, size(crop_10x,1)*px_10x);
fprintf('[CROP] GT nuclei in overlap: %d\n', numel(unique(crop_gt(crop_gt>0))));
%% -- Step 5: Intermediate diagnostic (pre-refinement) --------------------
plot_overlap_diagnostic(crop_10x, crop_gt, px_10x, params);
%% -- Step 5.5: Large-template NCC refinement -----------------------------
% Stage coordinates drift between sessions; the initial placement can be
% off by tens of um.  Instead of sliding small windows, we use the ENTIRE
% 10X segmentation mask as a single template and correlate against the
% padded 40X GT mask.  The multi-cell cluster provides many simultaneous
% correspondences which uniquely constrain the global shift.
seg_params = params;
seg_params.nucleus_channel = 1;
seg_params.metadata.px_size_x = px_10x;
fprintf('\n[ALIGN] Initial 10X segmentation for NCC template...\n');
[~, seg_initial] = segment_nuclei(crop_10x, seg_params);
fprintf('[ALIGN] Template: %d initial 10X nuclei\n', ...
        numel(unique(seg_initial(seg_initial>0))));
if ~isfield(params, 'align_ncc_search_um'); params.align_ncc_search_um = 100; end
search_px = round(params.align_ncc_search_um / px_10x);
[shift_r, shift_c, peak_ncc] = refine_alignment_ncc( ...
    seg_initial > 0, crop_gt > 0, search_px);
fprintf('[ALIGN] Shift: (dr=%+d, dc=%+d) px = (%+.1f, %+.1f) um  peak NCC=%.4f\n', ...
        shift_r, shift_c, shift_r*px_10x, shift_c*px_10x, peak_ncc);
if shift_r ~= 0 || shift_c ~= 0
    crop_gt = imtranslate(crop_gt, [shift_c, shift_r], ...
                          'FillValues', 0, 'Method', 'nearest');
    fprintf('[ALIGN] GT mask shifted.  Post-shift GT cells in view: %d\n', ...
            numel(unique(crop_gt(crop_gt>0))));
    plot_overlap_diagnostic(crop_10x, crop_gt, px_10x, params, 'refined');
end
%% -- Step 5.6: Manual GT alignment + scoring-ROI selection ----------------
% Since padding revealed the full 40X FOV, large black regions now appear
% where the 10X has no coverage.  Scoring over those regions would
% produce false misses.  The manual-alignment window now also carries a
% persistent, editable orange rectangle: the user aligns the GT overlay,
% drags/resizes the rectangle to bound the region where BOTH 10X has
% signal AND 40X has labelled cells, and confirms once.  Alignment and
% scoring-ROI are committed together, replacing the old Step 5.7.
%
%  Controls (inside the window):
%    - Drag inside image but OUTSIDE the orange rectangle: shift GT overlay
%    - Drag the orange rectangle / its corner handles: move / resize the ROI
%    - Arrow keys: 1 px GT nudge;  Shift+arrow: 10 px nudge
%    - Enter / Confirm: commit current shift + ROI
%    - Esc / Skip: zero shift, no ROI crop
%    - R: reset shift to zero (ROI unchanged)
manual_align_enabled = ~isfield(params, 'manual_align_enabled') || params.manual_align_enabled;
roi_rect_10x = [];
if manual_align_enabled
    fprintf('\n[ALIGN] Opening manual alignment + ROI selection window...\n');
    [dr_manual, dc_manual, roi_rect_10x] = ...
        interactive_manual_align(crop_10x, crop_gt, px_10x, params);
    if dr_manual ~= 0 || dc_manual ~= 0
        crop_gt = imtranslate(crop_gt, [dc_manual, dr_manual], ...
                              'FillValues', 0, 'Method', 'nearest');
        fprintf('[ALIGN] Manual shift applied: (dr=%+d, dc=%+d) px = (%+.1f, %+.1f) um\n', ...
                dr_manual, dc_manual, dr_manual*px_10x, dc_manual*px_10x);
        plot_overlap_diagnostic(crop_10x, crop_gt, px_10x, params, 'manual_refined');
    else
        fprintf('[ALIGN] No manual shift applied.\n');
    end
    if ~isempty(roi_rect_10x)
        fprintf('[ROI]  Scoring ROI set in alignment window: x=%d y=%d w=%d h=%d (%.0f x %.0f um)\n', ...
                roi_rect_10x(1), roi_rect_10x(2), roi_rect_10x(3), roi_rect_10x(4), ...
                roi_rect_10x(3)*px_10x, roi_rect_10x(4)*px_10x);
    else
        fprintf('[ROI]  No scoring ROI drawn in alignment window.\n');
    end
else
    fprintf('[ALIGN] Manual alignment disabled (manual_align_enabled=false).\n');
end
%% -- Step 5.7: Fallback ROI selection (used only if Step 5.6 didn't set one)
% The combined alignment + ROI window in Step 5.6 normally produces the
% scoring rectangle.  This fallback only runs when manual alignment was
% disabled (manual_align_enabled=false) or the user skipped ROI drawing
% there, so backwards compatibility is preserved.
roi_enabled = ~isfield(params, 'gt_roi_enabled') || params.gt_roi_enabled;
if roi_enabled && isempty(roi_rect_10x)
    fprintf('\n[ROI]  Opening interactive ROI selection on 10X overlap image...\n');
    fprintf('[ROI]  Draw a rectangle over the valid 10X/40X overlap zone.\n');
    fprintf('[ROI]  Double-click INSIDE the rectangle to confirm.\n');
    fprintf('[ROI]  Close the window without drawing to skip (use full crop).\n\n');
    % -- Build rich overlay for the selection figure -----------------------
    % Panel layout mirrors plot_overlap_diagnostic Panel 3 (the one shown
    % in the screenshot the user referenced): 10X grayscale + red GT edges.
    % We also add cyan GT label centroids so the boundary of coverage is
    % immediately visible, and a yellow dashed guide showing the current
    % crop extent.
    gt_bw_sel  = crop_gt > 0;
    gt_edge_sel = edge(gt_bw_sel, 'canny');
    roi_rgb    = repmat(mat2gray(crop_10x), [1 1 3]);
    roi_rgb(:,:,1) = min(1, roi_rgb(:,:,1) + gt_edge_sel);   % red GT edges
    % Compute GT label centroids for annotation
    gt_ids_sel   = unique(crop_gt(crop_gt > 0));
    gt_cents_sel = regionprops(crop_gt, 'Centroid');
    roi_fig = figure('Name', ...
        'Option 7 -- ROI Selection on 10X FOV (draw on 10X image, double-click to confirm)', ...
        'NumberTitle', 'off', 'Units', 'normalized', ...
        'Position', [0.04 0.04 0.72 0.88]);
    ax_sel = axes(roi_fig, 'Units', 'normalized', 'Position', [0.05 0.08 0.90 0.82]);
    imagesc(ax_sel, roi_rgb);
    axis(ax_sel, 'image'); axis(ax_sel, 'off');
    hold(ax_sel, 'on');
    % Mark each GT cell centroid with cyan circle + label ID
    for k = 1:numel(gt_ids_sel)
        if gt_ids_sel(k) <= numel(gt_cents_sel) && ~isempty(gt_cents_sel(gt_ids_sel(k)).Centroid)
            cx_gt = gt_cents_sel(gt_ids_sel(k)).Centroid(1);
            cy_gt = gt_cents_sel(gt_ids_sel(k)).Centroid(2);
            plot(ax_sel, cx_gt, cy_gt, 'co', ...
                 'MarkerSize', 9, 'LineWidth', 1.5, 'MarkerFaceColor', 'none');
            text(ax_sel, cx_gt + 8, cy_gt, num2str(gt_ids_sel(k)), ...
                 'Color', 'cyan', 'FontSize', 6, 'FontWeight', 'bold');
        end
    end
    % Highlight edges of the current crop extent with a yellow dashed border
    [Hcrop, Wcrop] = size(crop_10x);
    plot(ax_sel, [1 Wcrop Wcrop 1 1], [1 1 Hcrop Hcrop 1], ...
         'y--', 'LineWidth', 1.5);
    text(ax_sel, Wcrop * 0.02, Hcrop * 0.04, ...
         sprintf('Full crop: %d +/- %d px  (%.0f +/- %.0f um)', ...
                 Wcrop, Hcrop, Wcrop * px_10x, Hcrop * px_10x), ...
         'Color', 'yellow', 'FontSize', 8, 'BackgroundColor', [0 0 0 0.45]);
    title(ax_sel, ...
        {['10X Mosaic Overlap  |  NCC-refined  |  ' ...
          sprintf('%d GT nuclei (cyan circles + red boundaries)', numel(gt_ids_sel))], ...
         'Draw a rectangle to define the valid scoring zone, then DOUBLE-CLICK inside it', ...
         'Close this window to skip ROI and score the full crop'}, ...
        'FontSize', 10, 'FontWeight', 'bold');
    % -- Let user draw the rectangle ---------------------------------------
    try
        h_rect = drawrectangle(ax_sel, ...
                     'Color',          [1.0 0.55 0.0], ...
                     'LineWidth',       2.5, ...
                     'Label',          'Valid Scoring ROI -- double-click to confirm', ...
                     'LabelTextColor', [1 1 1], ...
                     'FaceAlpha',       0.07);
        % Block until user double-clicks inside the rectangle
        wait(h_rect);
        pos = h_rect.Position;   % [x, y, w, h] in axes/pixel units
        roi_rect_10x = round(pos);
        % Clamp to crop_10x bounds
        roi_rect_10x(1) = max(1, roi_rect_10x(1));
        roi_rect_10x(2) = max(1, roi_rect_10x(2));
        roi_rect_10x(3) = min(roi_rect_10x(3), Wcrop - roi_rect_10x(1) + 1);
        roi_rect_10x(4) = min(roi_rect_10x(4), Hcrop - roi_rect_10x(2) + 1);
        fprintf('[ROI]  Confirmed: x=%d  y=%d  w=%d  h=%d  (%.0f +/- %.0f um)\n', ...
                roi_rect_10x(1), roi_rect_10x(2), ...
                roi_rect_10x(3), roi_rect_10x(4), ...
                roi_rect_10x(3) * px_10x, roi_rect_10x(4) * px_10x);
        % Redraw confirmed box in solid orange over the dashed guide
        rectangle(ax_sel, 'Position', roi_rect_10x, ...
                  'EdgeColor', [1 0.55 0], 'LineWidth', 3.0);
        title(ax_sel, ...
            {sprintf('ROI confirmed: x=%d  y=%d  w=%d  h=%d  (%.0f +/- %.0f um)', ...
                     roi_rect_10x(1), roi_rect_10x(2), ...
                     roi_rect_10x(3), roi_rect_10x(4), ...
                     roi_rect_10x(3) * px_10x, roi_rect_10x(4) * px_10x), ...
             'Orange solid box = scoring zone  |  Segmentation + ROC run inside this region only'}, ...
            'FontSize', 10, 'FontWeight', 'bold');
        drawnow;
    catch ME_roi
        warning('[ROI]  drawrectangle cancelled or failed (%s). Using full crop.', ME_roi.message);
        roi_rect_10x = [];
    end
    % Save the ROI figure regardless of whether a rect was drawn
    if isfield(params, 'log_dir') && ~isempty(params.log_dir)
        roi_fig_path = fullfile(params.log_dir, ...
            sprintf('roi_selection_10x_%s.png', params.run_timestamp));
        exportgraphics(roi_fig, roi_fig_path, 'Resolution', 150);
        fprintf('[ROI]  Selection figure saved -> roi_selection_10x_%s.png\n', ...
                params.run_timestamp);
    end
    pause(1.5);
    close(roi_fig);
elseif ~isempty(roi_rect_10x)
    fprintf('[ROI]  Using ROI already set in Step 5.6; fallback selection skipped.\n');
else
    fprintf('[ROI]  ROI selection disabled (gt_roi_enabled=false). Scoring full crop.\n');
end
% -- Apply ROI crop to both crop_10x and crop_gt --------------------------
if ~isempty(roi_rect_10x)
    rx = roi_rect_10x(1);  ry = roi_rect_10x(2);
    rw = roi_rect_10x(3);  rh = roi_rect_10x(4);
    col_roi = rx : (rx + rw - 1);
    row_roi = ry : (ry + rh - 1);
    crop_10x = crop_10x(row_roi, col_roi);
    crop_gt  = crop_gt(row_roi,  col_roi);
    % Re-compact GT label IDs after crop (edge-clipped cells get dropped)
    gt_ids_after   = unique(crop_gt(crop_gt > 0));
    crop_gt_clean  = zeros(size(crop_gt), 'uint16');
    for k = 1:numel(gt_ids_after)
        crop_gt_clean(crop_gt == gt_ids_after(k)) = k;
    end
    crop_gt = crop_gt_clean;
    fprintf('[ROI]  After crop: %d +/- %d px  |  %d GT nuclei retained\n', ...
            rw, rh, numel(gt_ids_after));
    % Update crop_bbox to reflect the further sub-crop within the mosaic
    crop_bbox(1) = crop_bbox(1) + rx - 1;
    crop_bbox(2) = crop_bbox(2) + ry - 1;
    crop_bbox(3) = rw;
    crop_bbox(4) = rh;
    % Show a quick confirmation diagnostic of the final scoring region
    plot_overlap_diagnostic(crop_10x, crop_gt, px_10x, params, 'roi_cropped');
end
%% -- Step 6: Run segmentation methods on 10X overlap ---------------------
fprintf('\n[SEG]  Running segmentation on 10X overlap...\n');
% Method 1: Otsu baseline (L1) -- lower-tier, no if-then rules
fprintf('[SEG]  Method 1: Otsu baseline (L1)...\n');
[~, seg_matlab] = segment_nuclei_l1(crop_10x, seg_params);
fprintf('[SEG]    Detected %d connected components\n', max(seg_matlab(:)));
% Method 2: Rule-based pipeline + circularity filter (L2)
fprintf('[SEG]  Method 2: Rule-based + circularity filter (L2)...\n');
[cells_rule, seg_rule] = segment_nuclei(crop_10x, seg_params);
circ_thresh  = 0.60;
circ_vals    = [cells_rule.circularity];
pass         = circ_vals >= circ_thresh;
seg_L2       = zeros(size(seg_rule), 'uint16');
cids         = find(pass);
for k = 1:numel(cids); seg_L2(seg_rule == cids(k)) = k; end
fprintf('[SEG]    thresh=%.2f  %d/%d cells pass circularity\n', ...
        circ_thresh, sum(pass), numel(cells_rule));
% Method 3: StarDist
fprintf('[SEG]  Method 3: StarDist (L3)...\n');
seg_stardist = run_stardist_on_image(crop_10x);
fprintf('[SEG]    StarDist detected %d nuclei\n', max(seg_stardist(:)));
%% -- Step 7: Compute ROC vs 40X GT ---------------------------------------
fprintf('\n[ROC]  Computing IoU-based ROC vs 40X GT mask...\n');
methods = {'Otsu (L1)', 'Rule-based + Circ. (L2)', 'StarDist (L3)'};
segs    = {seg_matlab, seg_L2, seg_stardist};
colors  = {[0.3 0.6 1.0], [1.0 0.4 0.3], [0.2 0.8 0.3]};
roc     = cell(1,3);
for m = 1:3
    roc{m} = compute_roc_objectlevel(crop_gt, segs{m});
    % Real F1 at IoU=0.5 using ALL predictions (not the old score-gated one)
    P = roc{m}.tp_at_half / max(1, roc{m}.n_pred);
    R = roc{m}.tp_at_half / max(1, roc{m}.n_gt);
    roc{m}.f1_iou05 = 2*P*R / (P + R + eps);
    fprintf('[ROC]  %s: AP=%.3f  F1@IoU=0.5=%.3f  TP=%d/%d (rec=%.2f)  pred=%d (prec=%.2f)\n', ...
            methods{m}, roc{m}.ap, roc{m}.f1_iou05, ...
            roc{m}.tp_at_half, roc{m}.n_gt, R, roc{m}.n_pred, P);
end
%% -- Step 8: Display results ----------------------------------------------
plot_validation_results(crop_10x, crop_gt, segs, roc, methods, colors, px_10x, params, roi_rect_10x);
%% -- Step 9 (optional): Export centroid CSV + accuracy histogram ---------
% Mirrors Option 6 export. Gated by params.gt_export_centroids_csv.
export_on = isfield(params, 'gt_export_centroids_csv') && params.gt_export_centroids_csv;
acc_data  = [];
if export_on
    fprintf('\n[VAL]  Exporting centroid CSVs and accuracy histogram...\n');
    method_tags = {'L1', 'L2', 'L3'};
    acc_data    = struct();
    for m = 1:3
        seg      = segs{m};
        pred_ids = unique(seg(seg > 0));
        n_pred   = numel(pred_ids);
        iou_mat  = roc{m}.iou_matrix;
        gt_ids   = unique(crop_gt(crop_gt > 0));
        n_gt     = numel(gt_ids);
        % Greedy IoU matching at threshold 0.5
        matched_gt   = false(n_gt, 1);
        matched_pred = false(n_pred, 1);
        if n_pred > 0 && size(iou_mat,1) == n_gt && size(iou_mat,2) == n_pred
            [sorted_iou, sort_idx] = sort(iou_mat(:), 'descend');
            for si_v = 1:numel(sorted_iou)
                if sorted_iou(si_v) < 0.5; break; end
                [gi, pj] = ind2sub([n_gt, n_pred], sort_idx(si_v));
                if ~matched_gt(gi) && ~matched_pred(pj)
                    matched_gt(gi)   = true;
                    matched_pred(pj) = true;
                end
            end
        end
        % Per-prediction centroid table (in 10X mosaic crop coords)
        cell_id_out  = zeros(n_pred, 1);
        centroid_x_px= zeros(n_pred, 1);
        centroid_y_px= zeros(n_pred, 1);
        centroid_x_um= zeros(n_pred, 1);
        centroid_y_um= zeros(n_pred, 1);
        matched_gt_id= zeros(n_pred, 1);
        iou_score    = zeros(n_pred, 1);
        for pj = 1:n_pred
            mask_p = (seg == pred_ids(pj));
            props  = regionprops(mask_p, 'Centroid');
            if ~isempty(props)
                cx = props(1).Centroid(1);
                cy = props(1).Centroid(2);
            else
                cx = 0; cy = 0;
            end
            cell_id_out(pj)   = pred_ids(pj);
            centroid_x_px(pj) = cx;
            centroid_y_px(pj) = cy;
            centroid_x_um(pj) = cx * px_10x;
            centroid_y_um(pj) = cy * px_10x;
            if size(iou_mat,2) >= pj
                [best_iou, best_gi] = max(iou_mat(:, pj));
                iou_score(pj)       = best_iou;
                if matched_pred(pj)
                    matched_gt_id(pj) = gt_ids(best_gi);
                end
            end
        end
        T = table(cell_id_out, centroid_x_px, centroid_y_px, ...
                  centroid_x_um, centroid_y_um, ...
                  matched_gt_id, iou_score, ...
                  'VariableNames', {'cell_id', ...
                                    'centroid_x_px', 'centroid_y_px', ...
                                    'centroid_x_um', 'centroid_y_um', ...
                                    'matched_gt_id', 'iou_score'});
        csv_out = fullfile(params.log_dir, ...
            sprintf('centroids_10x_%s_%s.csv', method_tags{m}, params.run_timestamp));
        writetable(T, csv_out);
        fprintf('[VAL]  %s: centroid CSV saved (%d cells) -> %s\n', ...
                method_tags{m}, n_pred, csv_out);
        % Detection rate vs IoU threshold for histogram panel
        iou_bins = 0:0.05:1.0;
        det_pct  = zeros(size(iou_bins));
        if n_pred > 0 && size(iou_mat,2) == n_pred && n_gt > 0
            best_iou_per_gt = max(iou_mat, [], 2);
            for bi = 1:numel(iou_bins)
                det_pct(bi) = sum(best_iou_per_gt >= iou_bins(bi)) / n_gt * 100;
            end
        end
        acc_data(m).label    = methods{m};
        acc_data(m).tag      = method_tags{m};
        acc_data(m).iou_bins = iou_bins;
        acc_data(m).det_pct  = det_pct;
        acc_data(m).n_gt     = n_gt;
        acc_data(m).n_pred   = n_pred;
        acc_data(m).tp       = sum(matched_pred);
        acc_data(m).fp       = sum(~matched_pred);
        acc_data(m).fn       = sum(~matched_gt);
    end
    plot_accuracy_histogram(acc_data, roc, methods, colors, params);
end
%% -- Package output -------------------------------------------------------
results.mosaic_10x     = mosaic_10x;
results.gt_mask_canvas = gt_mask_canvas;
results.crop_10x       = crop_10x;
results.crop_gt        = crop_gt;
results.crop_bbox      = crop_bbox;
results.segs           = segs;
results.roc            = roc;
results.methods        = methods;
results.canvas_info    = canvas_info;
results.fov_40x_coords = fov_40x_coords;
results.px_10x         = px_10x;
results.px_40x         = px_40x;
results.roi_rect_10x   = roi_rect_10x;   % [] = full crop; [x y w h] = user-defined ROI
if export_on
    results.acc_data = acc_data;
end
end
%% ======================================================================
function f1 = get_f1_iou05(r)
%GET_F1_IOU05  Real F1 at IoU=0.5 computed from tp_at_half / n_pred / n_gt.
P  = r.tp_at_half / max(1, r.n_pred);
R  = r.tp_at_half / max(1, r.n_gt);
f1 = 2*P*R / (P + R + eps);
end
%% ======================================================================
function [cells, seg_mask, bw_nuclei] = segment_nuclei_l1(nucleus_img, params) %#ok<INUSD>
%SEGMENT_NUCLEI_L1  Otsu + connected-components baseline (Level 1).
fprintf('[SEG-L1] *** Otsu baseline IS being called (graythresh + bwlabel)\n');
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
%% ======================================================================
function [mosaic, canvas_info, px_10x] = build_10x_mosaic(params)
%BUILD_10X_MOSAIC  Stitch all 10X FOVs using stage coordinates.
data = bfopen(params.nd2_10x_path);
n_series = size(data, 1);
omeMeta  = data{1, 4};
try
    px_10x = double(omeMeta.getPixelsPhysicalSizeX(0).value());
catch
    px_10x = 0.655;
end
% Read stage positions for each series
stage_x = zeros(1, n_series);
stage_y = zeros(1, n_series);
for s = 1:n_series
    try
        stage_x(s) = double(omeMeta.getPlanePositionX(s-1, 0).value());
        stage_y(s) = double(omeMeta.getPlanePositionY(s-1, 0).value());
    catch
        warning('Could not read stage position for 10X series %d', s);
    end
end
% cameraTransformationMatrix: stage +X maps to image LEFT (from overlap_summary)
% So to convert stage coords to pixel coords, we invert the X direction
[H_tile, W_tile] = size(data{1,1}{params.channel_10x, 1});
% Convert stage positions to canvas pixel coords
% Each FOV centred at (stage_x_s, stage_y_s) in um
% Canvas origin: min stage_x, min stage_y corner
stage_x_img = -stage_x;   % invert due to camera transform
stage_y_img = -stage_y;   % invert Y (image convention vs stage convention)
stage_x_px = (stage_x_img - min(stage_x_img)) / px_10x;
stage_y_px = (stage_y_img - min(stage_y_img)) / px_10x;
% Canvas dimensions
canvas_W = round(max(stage_x_px) + W_tile);
canvas_H = round(max(stage_y_px) + H_tile);
mosaic   = zeros(canvas_H, canvas_W);
% Place each tile
for s = 1:n_series
    raw = double(data{s,1}{params.channel_10x, 1});
    raw = raw ./ (max(raw(:)) + eps);
    % Top-left corner of this tile in canvas
    tl_c = round(stage_x_px(s)) + 1;
    tl_r = round(stage_y_px(s)) + 1;
    r1 = tl_r; r2 = tl_r + H_tile - 1;
    c1 = tl_c; c2 = tl_c + W_tile - 1;
    r1 = max(1,r1); r2 = min(canvas_H, r2);
    c1 = max(1,c1); c2 = min(canvas_W, c2);
    tile_h = r2 - r1 + 1;
    tile_w = c2 - c1 + 1;
    mosaic(r1:r2, c1:c2) = max(mosaic(r1:r2, c1:c2), raw(1:tile_h, 1:tile_w));
end
canvas_info.stage_x_px = stage_x_px;
canvas_info.stage_y_px = stage_y_px;
canvas_info.stage_x_um = stage_x;
canvas_info.stage_y_um = stage_y;
canvas_info.tile_H     = H_tile;
canvas_info.tile_W     = W_tile;
canvas_info.canvas_H   = canvas_H;
canvas_info.canvas_W   = canvas_W;
canvas_info.x_origin_um = min(stage_x_img);
canvas_info.y_origin_um = min(stage_y_img);
fprintf('[10X]  %d FOVs placed on %dx%d canvas\n', n_series, canvas_W, canvas_H);
end
%% ======================================================================
function [mips, coords, px_40x] = load_40x_mips(params)
%LOAD_40X_MIPS  Load 40X FOVs, compute MIP per FOV, read stage coords.
data = bfopen(params.nd2_40x_path);
n_series_total = size(data, 1);
omeMeta = data{1, 4};
try
    px_40x = double(omeMeta.getPixelsPhysicalSizeX(0).value());
catch
    px_40x = 0.1625;
end
zpf     = params.z_planes_per_40x_fov;
n_fov   = floor(n_series_total / zpf);
[H40, W40] = size(data{1,1}{params.channel_40x, 1});
mips   = cell(1, n_fov);
coords = struct('fov', num2cell(1:n_fov), ...
                'stage_x_um', num2cell(zeros(1,n_fov)), ...
                'stage_y_um', num2cell(zeros(1,n_fov)));
for f = 1:n_fov
    s_start = (f-1)*zpf + 1;
    s_end   = s_start + zpf - 1;
    % MIP across Z stack for this FOV
    mip = zeros(H40, W40);
    for s = s_start:s_end
        plane = double(data{s,1}{params.channel_40x, 1});
        mip   = max(mip, plane);
    end
    mips{f} = mip ./ (max(mip(:)) + eps);
    % Stage coords of centre series
    mid_s = s_start + floor(zpf/2);
    try
        coords(f).stage_x_um = double(omeMeta.getPlanePositionX(mid_s-1, 0).value());
        coords(f).stage_y_um = double(omeMeta.getPlanePositionY(mid_s-1, 0).value());
    catch
    end
    fprintf('[40X]  FOV %d: stage = (%.1f, %.1f) um\n', f, ...
            coords(f).stage_x_um, coords(f).stage_y_um);
end
end
%% ======================================================================
function [gt_canvas, fov_rect_canvas, mosaic_out, canvas_info] = ...
    place_gt_in_canvas(params, fov_40x_coords, canvas_info, mosaic_in, px_10x, px_40x)
%PLACE_GT_IN_CANVAS  Place 40X GT mask into a 10X canvas that is padded,
%  if necessary, so the full 40X FOV fits without clipping.
%
%  Returns
%    gt_canvas        uint16 label image the size of the (possibly padded)
%                     10X canvas, with the full 40X GT mask pasted at its
%                     physical position -- no cells dropped at the edges.
%    fov_rect_canvas  [x y w h] in (padded) canvas coords: the full 40X
%                     FOV rectangle, unclipped.
%    mosaic_out       10X mosaic, zero-padded on whichever sides the 40X
%                     FOV extended beyond the original canvas.  Shares the
%                     coordinate frame of gt_canvas / fov_rect_canvas so
%                     crop_overlap can use them together.
%    canvas_info      updated with new canvas_H/W and origin_um offsets.
% Load GT mask (Frame 0 = nuclei)
gt_mask = uint16(imread(params.gt_mask_path, 1));
% Scale GT mask to 10X pixel size (nearest-neighbour to preserve labels)
scale_factor = px_40x / px_10x;
gt_scaled    = imresize(gt_mask, scale_factor, 'nearest');
[H_gs, W_gs] = size(gt_scaled);
n_cells_in_mask = numel(unique(gt_scaled(gt_scaled > 0)));
% Place at FOV1 location in canvas (apply same transform as 10X)
fov1_stage_x_img = -fov_40x_coords(1).stage_x_um;
fov1_stage_y_img = -fov_40x_coords(1).stage_y_um;
% Convert to canvas pixel coords using canvas origin
fov1_cx_px = (fov1_stage_x_img - canvas_info.x_origin_um) / px_10x;
fov1_cy_px = (fov1_stage_y_img - canvas_info.y_origin_um) / px_10x;
% Unclipped top-left / bottom-right of the 40X FOV in the original canvas
tl_c = round(fov1_cx_px - W_gs/2 + canvas_info.tile_W/2) + 1;
tl_r = round(fov1_cy_px - H_gs/2 + canvas_info.tile_H/2) + 1;
br_c = tl_c + W_gs - 1;
br_r = tl_r + H_gs - 1;
% How much padding is needed on each side to fit the full FOV?
pad_top    = max(0, 1 - tl_r);
pad_left   = max(0, 1 - tl_c);
pad_bottom = max(0, br_r - canvas_info.canvas_H);
pad_right  = max(0, br_c - canvas_info.canvas_W);
mosaic_out = mosaic_in;
if pad_top + pad_left + pad_bottom + pad_right > 0
    fprintf('[GT]   40X FOV extends beyond 10X canvas -- padding so full mask fits.\n');
    fprintf('[GT]     Pad (px): top=%d  left=%d  bottom=%d  right=%d\n', ...
            pad_top, pad_left, pad_bottom, pad_right);
    mosaic_out = padarray(mosaic_out, [pad_top, pad_left], 0, 'pre');
    mosaic_out = padarray(mosaic_out, [pad_bottom, pad_right], 0, 'post');
    canvas_info.canvas_H    = size(mosaic_out, 1);
    canvas_info.canvas_W    = size(mosaic_out, 2);
    canvas_info.x_origin_um = canvas_info.x_origin_um - pad_left * px_10x;
    canvas_info.y_origin_um = canvas_info.y_origin_um - pad_top  * px_10x;
    % Shift FOV placement into padded coords
    tl_c = tl_c + pad_left;   br_c = br_c + pad_left;
    tl_r = tl_r + pad_top;    br_r = br_r + pad_top;
end
% With padding applied, the full FOV now fits entirely -- no clipping
gt_canvas = zeros(canvas_info.canvas_H, canvas_info.canvas_W, 'uint16');
gt_canvas(tl_r:br_r, tl_c:br_c) = gt_scaled;
fov_rect_canvas = [tl_c, tl_r, W_gs, H_gs];
n_cells_placed = numel(unique(gt_canvas(gt_canvas > 0)));
fprintf('[GT]   GT mask (%dx%d, %d cells) placed at canvas (%d,%d)\n', ...
        W_gs, H_gs, n_cells_in_mask, tl_c, tl_r);
fprintf('[GT]   40X FOV rect in canvas: x=%d y=%d w=%d h=%d (%.0f x %.0f um)\n', ...
        fov_rect_canvas(1), fov_rect_canvas(2), ...
        fov_rect_canvas(3), fov_rect_canvas(4), ...
        fov_rect_canvas(3) * px_10x, fov_rect_canvas(4) * px_10x);
fprintf('[GT]   Cells actually placed: %d / %d\n', n_cells_placed, n_cells_in_mask);
end
%% ======================================================================
function [crop_img, crop_mask, bbox] = crop_overlap(mosaic, gt_canvas, fov_rect_canvas)
%CROP_OVERLAP  Extract the physical 40X FOV region from the 10X mosaic.
%
%  fov_rect_canvas (optional): [x y w h] full 40X FOV rectangle in canvas
%    pixel coords (from place_gt_in_canvas).  When provided, the crop is
%    the full 40X imaging area intersected with the 10X mosaic's imaged
%    region -- this keeps every 10X nucleus that lies under the 40X
%    objective in view, even those not annotated in the GT mask.
%
%  When fov_rect_canvas is omitted (legacy), the crop falls back to the
%  bounding box of GT-labelled cells.  This is kept so older callers don't
%  break, but new callers should always pass the rectangle.
if nargin >= 3 && ~isempty(fov_rect_canvas)
    c1 = fov_rect_canvas(1);
    r1 = fov_rect_canvas(2);
    c2 = c1 + fov_rect_canvas(3) - 1;
    r2 = r1 + fov_rect_canvas(4) - 1;
    fprintf('[CROP] Using 40X FOV rectangle: x=%d y=%d w=%d h=%d\n', ...
            c1, r1, c2 - c1 + 1, r2 - r1 + 1);
else
    [r_mask, c_mask] = find(gt_canvas > 0);
    if isempty(r_mask)
        error('No GT mask overlap with 10X canvas. Check stage coordinates.');
    end
    r1 = min(r_mask); r2 = max(r_mask);
    c1 = min(c_mask); c2 = max(c_mask);
    fprintf('[CROP] Legacy path: bbox of GT cells: x=%d y=%d w=%d h=%d\n', ...
            c1, r1, c2 - c1 + 1, r2 - r1 + 1);
end
% Clamp only to canvas bounds.  We deliberately do NOT clip to the
% mosaic's nonzero bounding box: the 40X FOV may extend slightly beyond
% the stitched 10X coverage, and any such uncovered area inside the crop
% simply shows black, letting the user see exactly where 10X data is
% missing.  The previous behaviour shrank the crop until every edge row
% had some 10X signal, which hid the rightmost / bottommost 10X cells
% the 40X covered but whose row/col of mosaic happened to be near-zero.
[Hm, Wm] = size(mosaic);
r1_raw = r1;  r2_raw = r2;  c1_raw = c1;  c2_raw = c2;
r1 = max(1, r1);   r2 = min(Hm, r2);
c1 = max(1, c1);   c2 = min(Wm, c2);
if r1 ~= r1_raw || r2 ~= r2_raw || c1 ~= c1_raw || c2 ~= c2_raw
    fprintf(['[CROP] FOV rect clipped to canvas bounds: ' ...
             '(%d,%d)-(%d,%d)  ->  (%d,%d)-(%d,%d)\n'], ...
             c1_raw, r1_raw, c2_raw, r2_raw, c1, r1, c2, r2);
end
if r2 < r1 || c2 < c1
    error('40X FOV rectangle does not intersect the 10X canvas.');
end
crop_img  = mosaic(r1:r2, c1:c2);
crop_mask = gt_canvas(r1:r2, c1:c2);
bbox      = [c1, r1, c2-c1+1, r2-r1+1];
fprintf('[CROP] Final crop size: %d x %d px\n', c2 - c1 + 1, r2 - r1 + 1);
end
%% ======================================================================
function seg = run_stardist_on_image(img)
%RUN_STARDIST_ON_IMAGE  Run StarDist via bridge on a 2D image.
seg = zeros(size(img), 'uint16');
bridge_path = fullfile(fileparts(mfilename('fullpath')), 'stardist_bridge.py');
if ~isfile(bridge_path); warning('Bridge not found'); return; end
try
    spec   = py.importlib.util.spec_from_file_location('stardist_bridge', bridge_path);
    bridge = py.importlib.util.module_from_spec(spec);
    spec.loader.exec_module(bridge);
    load_fn    = py.getattr(bridge, 'load_model');
    model      = feval(load_fn, '2D_versatile_fluo');
    predict_fn = py.getattr(bridge, 'predict_full');
    np         = py.importlib.import_module('numpy');
    img_f32    = single(img) ./ (max(single(img(:))) + eps);
    [H, W]     = size(img_f32);
    py_flat    = np.array(img_f32(:));
    result     = feval(predict_fn, model, py_flat, int32(H), int32(W));
    labels_list = result{'labels_flat'};
    labels_flat = double(py.array.array('l', labels_list));
    seg = uint16(reshape(labels_flat, H, W));
catch ME
    warning('StarDist failed: %s', ME.message);
end
end
%% ======================================================================
function roc = compute_roc_objectlevel(gt_mask, pred_mask)
%COMPUTE_ROC_OBJECTLEVEL  Score-threshold ROC using area as proxy score.
iou_fixed = 0.5;
gt_ids   = unique(gt_mask(gt_mask > 0));
pred_ids = unique(pred_mask(pred_mask > 0));
n_gt     = numel(gt_ids);
n_pred   = numel(pred_ids);
iou_matrix  = zeros(n_gt, n_pred);
pred_scores = zeros(1, n_pred);
for pi = 1:n_pred
    pm = (pred_mask == pred_ids(pi));
    pred_scores(pi) = sum(pm(:));
    for gi = 1:n_gt
        gm = (gt_mask == gt_ids(gi));
        inter = sum(gm(:) & pm(:));
        if inter > 0
            uni = sum(gm(:) | pm(:));
            iou_matrix(gi, pi) = inter / uni;
        end
    end
end
if max(pred_scores) > 0; pred_scores = pred_scores ./ max(pred_scores); end
score_thr = [1.0, sort(unique(pred_scores),'descend'), -eps];
n_t = numel(score_thr);
tpr = zeros(1,n_t); fpr = zeros(1,n_t); prec = zeros(1,n_t);
for ti = 1:n_t
    keep = pred_scores >= score_thr(ti);
    keep_idx = find(keep);
    n_keep   = numel(keep_idx);
    m_gt   = false(n_gt,1);
    m_pred = false(n_keep,1);
    if n_keep > 0
        sub = iou_matrix(:, keep_idx);
        [si, idx] = sort(sub(:), 'descend');
        for k = 1:numel(si)
            if si(k) < iou_fixed; break; end
            [gi, pi] = ind2sub([n_gt, n_keep], idx(k));
            if ~m_gt(gi) && ~m_pred(pi); m_gt(gi)=true; m_pred(pi)=true; end
        end
    end
    TP = sum(m_gt); FP = n_keep - sum(m_pred);
    tpr(ti) = TP / (n_gt + eps);
    fpr(ti) = FP / (FP + n_gt + eps);
    if TP+FP > 0; prec(ti) = TP/(TP+FP); end
end
[fs, idx] = sort(fpr);
ts = tpr(idx);
auc = trapz(fs, ts);
[~, mid] = min(abs(score_thr - 0.5));
f1 = 2*prec(mid)*tpr(mid) / (prec(mid) + tpr(mid) + eps);
roc.thresholds = score_thr;
roc.tpr = tpr; roc.fpr = fpr; roc.precision = prec;
roc.auc = auc; roc.f1_at_half = f1;
roc.n_gt = n_gt; roc.n_pred = n_pred;
roc.iou_matrix = iou_matrix;
% -- PR curve over IoU thresholds (COCO-style) ----------------------------
% Sweep IoU threshold from strict (0.95) down to lenient (0.05). At each
% threshold, greedily match predictions to GT (1:1, highest IoU first).
% Produces a proper curve for object detection where no calibrated
% confidence score exists on the prediction side.
iou_grid    = 0.05:0.05:0.95;
n_iou       = numel(iou_grid);
pr_prec     = zeros(1, n_iou);
pr_rec      = zeros(1, n_iou);
tp_at_half  = 0;
for ii = 1:n_iou
    thr = iou_grid(ii);
    m_gt   = false(n_gt,1);
    m_pred = false(n_pred,1);
    [si, idx] = sort(iou_matrix(:), 'descend');
    for k = 1:numel(si)
        if si(k) < thr; break; end
        [gi, pi] = ind2sub([n_gt, n_pred], idx(k));
        if ~m_gt(gi) && ~m_pred(pi)
            m_gt(gi) = true; m_pred(pi) = true;
        end
    end
    TP = sum(m_gt);
    if n_pred > 0; pr_prec(ii) = TP / n_pred; else; pr_prec(ii) = 1; end
    if n_gt   > 0; pr_rec(ii)  = TP / n_gt;   else; pr_rec(ii)  = 0; end
    if abs(thr - 0.5) < 1e-9; tp_at_half = TP; end
end
% Average precision = area under PR curve (sorted by ascending recall)
[rec_s, si] = sort(pr_rec);
prec_s      = pr_prec(si);
ap          = trapz([0 rec_s 1], [prec_s(1) prec_s prec_s(end)]);
roc.iou_grid     = iou_grid;
roc.pr_precision = pr_prec;
roc.pr_recall    = pr_rec;
roc.ap           = ap;
roc.tp_at_half   = tp_at_half;
end
%% ======================================================================
function [shift_r, shift_c, peak_ncc] = refine_alignment_ncc(seg_bin, gt_bin, search_px)
%REFINE_ALIGNMENT_NCC  Find integer shift (shift_r, shift_c) to apply to
%  gt_bin so that it best aligns with seg_bin.
%
%  Uses the WHOLE 10X segmentation mask as a single large template and
%  correlates against a padded GT mask (normxcorr2).  This is the
%  "large-window / whole-cluster" approach the user requested: one big
%  template with many simultaneous cell correspondences constrains the
%  global shift far more robustly than sliding small patches, which can
%  hunt for local optima in noisy regions.
[H, W] = size(seg_bin);
gt_pad = padarray(double(gt_bin), [search_px, search_px], 0, 'both');
C = normxcorr2(double(seg_bin), gt_pad);
% Zero-shift corresponds to template's top-left at (search_px+1, search_px+1)
% in gt_pad, i.e. C row,col = (search_px + H, search_px + W).
zr = search_px + H;
zc = search_px + W;
% Restrict search to +/- search_px around zero-shift (ignore wrap-around peaks).
mask = false(size(C));
r_lo = max(1, zr - search_px); r_hi = min(size(C,1), zr + search_px);
c_lo = max(1, zc - search_px); c_hi = min(size(C,2), zc + search_px);
mask(r_lo:r_hi, c_lo:c_hi) = true;
Cm = C; Cm(~mask) = -inf;
[peak_ncc, idx] = max(Cm(:));
[pr, pc] = ind2sub(size(C), idx);
% Shift to APPLY to gt_bin so it aligns with seg_bin:
% If NCC peaks when seg_bin's top-left overlays gt's position (1+dr_peak, 1+dc_peak),
% then gt's cells are displaced by +dr_peak rows / +dc_peak cols relative to seg.
% To re-align, shift gt by (-dr_peak, -dc_peak).
shift_r = zr - pr;
shift_c = zc - pc;
end
%% ======================================================================
function plot_overlap_diagnostic(crop_10x, crop_gt, px_10x, params, tag)
%PLOT_OVERLAP_DIAGNOSTIC  Intermediate check: does 10X crop align with GT?
if nargin < 5 || isempty(tag); tag = 'initial'; end
% Suppress diagnostic plot unless explicitly enabled. Gated so all four
% callers (pre-NCC, refined, manual_refined, roi_cropped) short-circuit.
if ~isfield(params, 'show_alignment_diagnostics') || ~params.show_alignment_diagnostics
    fprintf('[DIAG] Alignment diagnostic (%s) skipped -- set params.show_alignment_diagnostics=true to show.\n', tag);
    return;
end
fig = figure('Name',['Overlap Diagnostic - ' tag],'NumberTitle','off', ...
             'Units','normalized','Position',[0.05 0.1 0.9 0.75]);
sgtitle(sprintf('10X Crop vs 40X GT Mask Alignment Check (%s)', tag), ...
        'FontSize',13,'FontWeight','bold');
ax1 = subplot(1,3,1);
imagesc(ax1, crop_10x); colormap(ax1,'gray'); axis(ax1,'image'); axis(ax1,'off');
title(ax1, sprintf('10X crop (%dx%d px, %.0fx%.0f um)', ...
      size(crop_10x,2), size(crop_10x,1), ...
      size(crop_10x,2)*px_10x, size(crop_10x,1)*px_10x), 'FontSize',10);
ax2 = subplot(1,3,2);
imagesc(ax2, crop_gt); colormap(ax2,'jet'); axis(ax2,'image'); axis(ax2,'off');
title(ax2, sprintf('40X GT mask in 10X scale (%d nuclei)', ...
      numel(unique(crop_gt(crop_gt>0)))),'FontSize',10);
ax3 = subplot(1,3,3);
rgb = repmat(mat2gray(crop_10x), [1 1 3]);
gt_edge = edge(crop_gt > 0, 'canny');
rgb(:,:,1) = min(1, rgb(:,:,1) + gt_edge);
imagesc(ax3, rgb); axis(ax3,'image'); axis(ax3,'off');
title(ax3, 'Overlay: 10X grayscale + GT boundaries (red)','FontSize',10);
drawnow;
if isfield(params,'log_dir') && ~isempty(params.log_dir)
    exportgraphics(fig, fullfile(params.log_dir, ...
        sprintf('overlap_diagnostic_%s_%s.png', tag, params.run_timestamp)), ...
        'Resolution', 150);
end
end
%% ======================================================================
function plot_validation_results(crop_10x, crop_gt, segs, roc, methods, colors, px_10x, params, roi_rect_10x)
%PLOT_VALIDATION_RESULTS  Segmentation overlays + ROC + morphology bars.
if nargin < 9; roi_rect_10x = []; end
roi_str = ' (full crop)';
if ~isempty(roi_rect_10x)
    roi_str = sprintf(' (ROI: x=%d y=%d w=%d h=%d  |  %.0f?%.0f um)', ...
                      roi_rect_10x(1), roi_rect_10x(2), ...
                      roi_rect_10x(3), roi_rect_10x(4), ...
                      roi_rect_10x(3)*px_10x, roi_rect_10x(4)*px_10x);
end
fig = figure('Name','10X vs 40X GT Validation','NumberTitle','off', ...
             'Units','normalized','Position',[0.02 0.02 0.96 0.92]);
sgtitle(['Segmentation Validation  10X input vs 40X GT mask' roi_str], ...
        'FontSize', 12, 'FontWeight', 'bold');
method_short = {'Otsu (L1)', 'Rule-based (L2)', 'StarDist (L3)'};
% Panels 1-3: overlays
for m = 1:3
    ax = subplot(2,3,m);
    pred_edge = edge(segs{m} > 0, 'canny');
    gt_edge   = edge(crop_gt > 0, 'canny');
    rgb = repmat(mat2gray(crop_10x), [1 1 3]);
    rgb(:,:,1) = rgb(:,:,1) + pred_edge * colors{m}(1) * 0.8;
    rgb(:,:,2) = rgb(:,:,2) + pred_edge * colors{m}(2) * 0.8;
    rgb(:,:,3) = rgb(:,:,3) + pred_edge * colors{m}(3) * 0.8;
    rgb(:,:,:) = rgb(:,:,:) + repmat(gt_edge * 0.5, [1 1 3]);
    imagesc(ax, min(1, rgb)); axis(ax,'image'); axis(ax,'off');
    f1_iou05 = get_f1_iou05(roc{m});
    title(ax, sprintf('%s\nAP=%.3f  F1@IoU=0.5=%.3f  matched=%d/%d', ...
          method_short{m}, roc{m}.ap, f1_iou05, ...
          roc{m}.tp_at_half, roc{m}.n_gt), ...
          'FontSize', 9);
end
% Panel 4: Precision-Recall curves over IoU threshold sweep (COCO-style)
ax4 = subplot(2,3,4);
hold(ax4,'on');
for m = 1:3
    % Sort by recall ascending so the line draws left-to-right
    [rec_s, si] = sort(roc{m}.pr_recall);
    prec_s      = roc{m}.pr_precision(si);
    plot(ax4, rec_s, prec_s, '-o', 'Color', colors{m}, 'LineWidth', 2, ...
         'MarkerSize', 4, 'MarkerFaceColor', colors{m}, ...
         'DisplayName', sprintf('%s (AP=%.3f)', method_short{m}, roc{m}.ap));
    % Mark the IoU=0.5 operating point
    [~, i05] = min(abs(roc{m}.iou_grid - 0.5));
    plot(ax4, roc{m}.pr_recall(i05), roc{m}.pr_precision(i05), ...
         'p', 'MarkerFaceColor', colors{m}, 'MarkerEdgeColor', 'k', ...
         'MarkerSize', 10, 'HandleVisibility', 'off');
end
xlabel(ax4, 'Recall'); ylabel(ax4, 'Precision');
title(ax4, 'Precision-Recall (IoU sweep 0.05-0.95)', 'FontSize', 10);
legend(ax4, 'Location', 'southwest','FontSize',8);
xlim(ax4,[0 1]); ylim(ax4,[0 1.05]); grid(ax4,'on'); axis(ax4,'square');
text(ax4, 0.02, 0.08, 'star = IoU 0.5', 'FontSize', 7, 'Color', [0.3 0.3 0.3]);
% Panel 5: Detection summary -- expressed as % of GT cell count
ax5 = subplot(2,3,5);
bar_data = zeros(3, 3);
for m = 1:3
    n_gt_m = max(1, roc{m}.n_gt);
    bar_data(m,1) = 100;                                  % GT = 100% baseline
    bar_data(m,2) = 100 * roc{m}.n_pred    / n_gt_m;      % predicted as % of GT
    bar_data(m,3) = 100 * roc{m}.tp_at_half / n_gt_m;     % matched @ IoU=0.5 (recall %)
end
b = bar(ax5, bar_data');
for m = 1:3; b(m).FaceColor = colors{m}; end
set(ax5, 'XTickLabel', {'GT nuclei','Predicted','Matched @ IoU=0.5'}, 'XTick', 1:3);
% Annotate each bar with its percentage
for m = 1:3
    for g = 1:3
        text(ax5, g + b(m).XOffset, bar_data(m,g) + 2, ...
             sprintf('%.0f%%', bar_data(m,g)), ...
             'HorizontalAlignment','center', 'FontSize', 7, 'Color', [0.2 0.2 0.2]);
    end
end
ylabel(ax5, '% of GT cells');
ylim(ax5, [0 max(110, 1.15*max(bar_data(:)))]);
legend(ax5, method_short, 'Location', 'best', 'FontSize', 8);
title(ax5, 'Detection Rate (% of GT)', 'FontSize', 10);
grid(ax5,'on');
% Panel 6: Crop + GT overlay reference
ax6 = subplot(2,3,6);
rgb = repmat(mat2gray(crop_10x), [1 1 3]);
gt_edge = edge(crop_gt > 0, 'canny');
rgb(:,:,1) = min(1, rgb(:,:,1) + gt_edge);
imagesc(ax6, rgb); axis(ax6,'image'); axis(ax6,'off');
title(ax6, '10X + GT boundaries', 'FontSize', 10);
drawnow;
if isfield(params,'log_dir') && ~isempty(params.log_dir)
    exportgraphics(fig, fullfile(params.log_dir, ...
        sprintf('validation_10x_40x_%s.png', params.run_timestamp)), ...
        'Resolution', 150);
end
end
%% ======================================================================
function plot_accuracy_histogram(acc_data, roc, ~, colors, params)
%PLOT_ACCURACY_HISTOGRAM  Per-method detection-rate histogram vs IoU
%  threshold + IoU=0.5 paired bar + TP/FP/FN bar. Mirrors Option 6.
fig = figure('Name', 'Segmentation Accuracy vs Ground Truth (10X)', ...
             'NumberTitle', 'off', 'Units', 'normalized', ...
             'Position', [0.04 0.08 0.92 0.82]);
sgtitle({'Segmentation Accuracy vs Ground Truth (10X -> 40X)', ...
         'Y-axis: % of GT cells detected at IoU \geq threshold  |  X-axis: IoU threshold'}, ...
        'FontSize', 13, 'FontWeight', 'bold');
method_short = {'Otsu (L1)', 'Rule-based (L2)', 'StarDist (L3)'};
pair_titles  = {'L1 / GT', 'L2 / GT', 'L3 / GT'};
iou_display  = 0:0.1:1.0;
% Panels 1-3: per-method bars
for m = 1:3
    ax = subplot(2, 4, m);
    det_pct_disp = interp1(acc_data(m).iou_bins, acc_data(m).det_pct, ...
                           iou_display, 'linear', 0);
    det_pct_disp = max(0, det_pct_disp);
    bar(ax, iou_display, det_pct_disp, 0.7, ...
        'FaceColor', colors{m}, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    hold(ax, 'on');
    yline(ax, 100, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.8, ...
          'Label', 'GT (100%)', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    det_at_half = interp1(acc_data(m).iou_bins, acc_data(m).det_pct, 0.5, 'linear', 0);
    plot(ax, 0.5, det_at_half, 'kv', 'MarkerSize', 9, ...
         'MarkerFaceColor', colors{m}, 'LineWidth', 1.5);
    text(ax, 0.52, det_at_half + 3, sprintf('%.1f%%', det_at_half), ...
         'FontSize', 8, 'FontWeight', 'bold', 'Color', colors{m} * 0.7);
    tp = acc_data(m).tp; fp = acc_data(m).fp; fn = acc_data(m).fn;
    prec = tp / (tp + fp + eps) * 100;
    rec  = tp / (tp + fn + eps) * 100;
    f1   = 2 * prec * rec / (prec + rec + eps);
    ann_str = sprintf('TP=%d  FP=%d  FN=%d\nPrec=%.1f%%  Rec=%.1f%%\nF1=%.1f%%  AP=%.3f', ...
                      tp, fp, fn, prec, rec, f1, roc{m}.ap);
    text(ax, 0.55, 60, ann_str, 'FontSize', 7.5, 'Color', [0.2 0.2 0.2], ...
         'BackgroundColor', [0.97 0.97 0.97], 'EdgeColor', colors{m}, ...
         'LineWidth', 0.8, 'VerticalAlignment', 'bottom');
    xlabel(ax, 'IoU Threshold');
    ylabel(ax, 'Detection Rate (% of GT)');
    title(ax, {pair_titles{m}, method_short{m}}, 'FontSize', 10, 'FontWeight', 'bold');
    xlim(ax, [-0.05 1.05]); ylim(ax, [0 115]);
    set(ax, 'XTick', iou_display, 'XTickLabelRotation', 45);
    grid(ax, 'on');
end
% Panel 4: overlay of all three methods
ax4 = subplot(2, 4, 4);
hold(ax4, 'on');
leg_h = gobjects(3, 1);
for m = 1:3
    leg_h(m) = plot(ax4, acc_data(m).iou_bins, acc_data(m).det_pct, '-o', ...
                    'Color', colors{m}, 'LineWidth', 2.2, 'MarkerSize', 5, ...
                    'DisplayName', method_short{m});
end
yline(ax4, 100, 'k--', 'LineWidth', 1.5, 'Label', 'GT (100%)', ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 8);
xline(ax4, 0.5, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
      'Label', 'IoU=0.5', 'LabelVerticalAlignment', 'bottom', 'FontSize', 8);
xlabel(ax4, 'IoU Threshold'); ylabel(ax4, 'Detection Rate (% of GT)');
legend(ax4, leg_h, 'Location', 'southwest', 'FontSize', 8);
title(ax4, {'All Methods -- Overlay', 'detection rate vs IoU threshold'}, ...
      'FontSize', 10, 'FontWeight', 'bold');
xlim(ax4, [0 1]); ylim(ax4, [0 115]); grid(ax4, 'on');
% Panel 5-6: Paired accuracy at IoU=0.5
ax5 = subplot(2, 4, [5 6]);
bar_vals = zeros(3, 1);
for m = 1:3
    bar_vals(m) = interp1(acc_data(m).iou_bins, acc_data(m).det_pct, 0.5, 'linear', 0);
end
b = bar(ax5, 1:3, bar_vals, 0.55, 'EdgeColor', 'none');
b.FaceColor = 'flat';
for m = 1:3
    b.CData(m,:) = colors{m};
end
hold(ax5, 'on');
yline(ax5, 100, 'k--', 'LineWidth', 1.8, 'Label', 'GT 100%', ...
      'LabelHorizontalAlignment', 'right', 'FontSize', 9);
for m = 1:3
    text(ax5, m, bar_vals(m) + 2.5, sprintf('%.1f%%', bar_vals(m)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, ...
         'FontWeight', 'bold', 'Color', colors{m} * 0.75);
end
set(ax5, 'XTick', 1:3, 'XTickLabel', {'L1 / GT', 'L2 / GT', 'L3 / GT'}, 'FontSize', 10);
ylabel(ax5, 'Detection Rate at IoU \geq 0.5 (%)');
title(ax5, 'Paired Accuracy at IoU = 0.5  (L1/GT  +/-  L2/GT  +/-  L3/GT)', ...
      'FontSize', 11, 'FontWeight', 'bold');
ylim(ax5, [0 115]); grid(ax5, 'on');
% Panel 7-8: TP/FP/FN per method
ax6 = subplot(2, 4, [7 8]);
tp_fp_fn = zeros(3, 3);
for m = 1:3
    tp_fp_fn(m,1) = acc_data(m).tp;
    tp_fp_fn(m,2) = acc_data(m).fp;
    tp_fp_fn(m,3) = acc_data(m).fn;
end
b2 = bar(ax6, tp_fp_fn, 'grouped');
b2(1).FaceColor = [0.25 0.75 0.25];
b2(2).FaceColor = [1.00 0.35 0.25];
b2(3).FaceColor = [0.55 0.55 0.55];
for m = 1:3
    for k = 1:3
        text(ax6, b2(k).XEndPoints(m), b2(k).YEndPoints(m) + 0.3, ...
             num2str(tp_fp_fn(m,k)), 'HorizontalAlignment', 'center', ...
             'FontSize', 8, 'FontWeight', 'bold');
    end
end
set(ax6, 'XTick', 1:3, 'XTickLabel', {'L1', 'L2', 'L3'}, 'FontSize', 10);
ylabel(ax6, 'Cell Count');
legend(ax6, {'TP', 'FP', 'FN'}, 'Location', 'best', 'FontSize', 9);
title(ax6, 'TP / FP / FN per Method  (IoU \geq 0.5)', ...
      'FontSize', 11, 'FontWeight', 'bold');
grid(ax6, 'on');
drawnow;
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_path = fullfile(params.log_dir, ...
        sprintf('accuracy_histogram_10x_%s.png', params.run_timestamp));
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Accuracy histogram saved -> accuracy_histogram_10x_%s.png\n', ...
            params.run_timestamp);
end
end
%% ======================================================================
function [dr, dc, roi_rect] = interactive_manual_align(crop_10x, crop_gt, px_10x, params)
%INTERACTIVE_MANUAL_ALIGN  Align 40X GT overlay AND pick scoring ROI in one window.
%
%  The 10X grayscale image is the static canvas.  On top sit two
%  overlays: (1) the draggable 40X GT (red boundaries + cyan centroids)
%  for alignment, and (2) a persistent orange rectangle for the scoring
%  ROI.  Returns the integer shift (dr, dc) that crop_gt should be
%  translated by, plus roi_rect = [x y w h] in crop_10x coords bounding
%  the scoring zone (empty if the user didn't want one).
%
%  Controls:
%    - Drag inside image OUTSIDE the orange rectangle:  shift GT overlay
%    - Drag the orange rectangle / its corner handles:  move / resize ROI
%    - Arrow keys: 1 px GT shift; Shift+arrow: 10 px GT shift
%    - Enter / Confirm button: commit current shift + ROI
%    - Esc / Skip button:      zero shift, no ROI crop
%    - R / Reset button:       zero the shift, keep window open (ROI unchanged)
%    - Close window (X):       treat as confirm with current state
dr = 0;
dc = 0;
roi_rect = [];
[H, W]   = size(crop_10x);
base     = mat2gray(crop_10x);
gt_edge0 = edge(crop_gt > 0, 'canny');
gt_stats = regionprops(crop_gt, 'Centroid');
gt_ids   = unique(crop_gt(crop_gt > 0));
% Auto-suggest initial ROI: intersection of 10X-imaged bbox and GT-cells bbox.
% This lands the rectangle on the region that's most likely to yield a valid
% scoring zone; the user just nudges the corners to refine.
[r_img_i, c_img_i] = find(crop_10x > 0);
[r_gt_i,  c_gt_i]  = find(crop_gt > 0);
init_roi = [round(W*0.15), round(H*0.15), round(W*0.70), round(H*0.70)];
if ~isempty(r_img_i) && ~isempty(r_gt_i)
    rmi = max(min(r_img_i), min(r_gt_i));
    rma = min(max(r_img_i), max(r_gt_i));
    cmi = max(min(c_img_i), min(c_gt_i));
    cma = min(max(c_img_i), max(c_gt_i));
    if rma > rmi + 20 && cma > cmi + 20
        init_roi = [cmi, rmi, cma - cmi + 1, rma - rmi + 1];
    end
end
% Clamp to image bounds
init_roi(1) = max(1, init_roi(1));
init_roi(2) = max(1, init_roi(2));
init_roi(3) = min(W - init_roi(1) + 1, max(10, init_roi(3)));
init_roi(4) = min(H - init_roi(2) + 1, max(10, init_roi(4)));
% Collect valid centroid positions (label gaps possible pre-ROI-crop)
orig_cx = zeros(1, numel(gt_ids));
orig_cy = zeros(1, numel(gt_ids));
orig_id = zeros(1, numel(gt_ids));
n_cent  = 0;
for k = 1:numel(gt_ids)
    gid = gt_ids(k);
    if gid <= numel(gt_stats) && ~isempty(gt_stats(gid).Centroid)
        n_cent = n_cent + 1;
        orig_cx(n_cent) = gt_stats(gid).Centroid(1);
        orig_cy(n_cent) = gt_stats(gid).Centroid(2);
        orig_id(n_cent) = gid;
    end
end
orig_cx = orig_cx(1:n_cent);
orig_cy = orig_cy(1:n_cent);
orig_id = orig_id(1:n_cent);
% Drag state
is_dragging   = false;
drag_start_xy = [0 0];
drag_start_dr = 0;
drag_start_dc = 0;
% NOTE: CloseRequestFcn uses the src argument (@(src,~) delete(src)) rather
% than capturing `fig` from the workspace.  If we wrote @(~,~) delete(fig),
% the anonymous function is created BEFORE fig is assigned -- the capture is
% unresolved -- and any later close/close-all (including clear;close all;
% at the top of run_benchmark.m) throws "Unrecognized function or
% variable 'fig'" when it tries to fire the callback on a stale figure.
fig = figure('Name', ...
    'Option 7 -- Manual GT Alignment (drag image or arrow keys to move 40X overlay over 10X)', ...
    'NumberTitle', 'off', 'Units', 'normalized', ...
    'Position', [0.04 0.04 0.72 0.88], ...
    'KeyPressFcn', @on_key, ...
    'WindowButtonDownFcn',   @on_button_down, ...
    'WindowButtonUpFcn',     @on_button_up, ...
    'WindowButtonMotionFcn', @on_button_motion, ...
    'CloseRequestFcn',       @(src,~) delete(src));
ax = axes(fig, 'Units', 'normalized', 'Position', [0.05 0.15 0.90 0.78]);
% Static 10X grayscale background
h_bg = imagesc(ax, base);
colormap(ax, 'gray');
axis(ax, 'image'); axis(ax, 'off');
set(ax, 'XLim', [0.5, W+0.5], 'YLim', [0.5, H+0.5], ...
        'XLimMode', 'manual', 'YLimMode', 'manual');
hold(ax, 'on');
set(h_bg, 'HitTest', 'off', 'PickableParts', 'none');
% Red-edge overlay as an RGB image with AlphaData = edge mask.
% Translating the overlay = just updating XData/YData of this image.
red_rgb = cat(3, ones(H, W), zeros(H, W), zeros(H, W));
h_edge = image(ax, red_rgb, ...
               'AlphaData', double(gt_edge0) * 0.9, ...
               'XData', [1, W], 'YData', [1, H]);
set(h_edge, 'HitTest', 'off', 'PickableParts', 'none');
% Cyan centroid markers + label IDs (all shifted together)
h_cent = plot(ax, orig_cx, orig_cy, 'co', ...
              'MarkerSize', 9, 'LineWidth', 1.5, ...
              'HitTest', 'off', 'PickableParts', 'none');
h_txt = gobjects(1, n_cent);
for k = 1:n_cent
    h_txt(k) = text(ax, orig_cx(k) + 8, orig_cy(k), num2str(orig_id(k)), ...
                    'Color', 'cyan', 'FontSize', 6, 'FontWeight', 'bold', ...
                    'HitTest', 'off', 'PickableParts', 'none');
end
% Persistent editable scoring-ROI rectangle.  Its own handles handle the
% move/resize interaction within its footprint; clicks OUTSIDE it fall
% through to the figure's WindowButtonDownFcn (alignment drag).
h_roi = images.roi.Rectangle(ax, ...
    'Position',         init_roi, ...
    'Color',            [1.0 0.55 0.0], ...
    'StripeColor',      'w', ...
    'LineWidth',        2.5, ...
    'FaceAlpha',        0.08, ...
    'Label',            'Scoring ROI -- drag to move, corners to resize', ...
    'LabelTextColor',   [1 1 1], ...
    'LabelAlpha',       0.4);
h_title = title(ax, '', 'FontSize', 10, 'FontWeight', 'bold');
% Control buttons -- also hook KeyPressFcn so arrow keys work after a click
% (focus may move to a button and would otherwise swallow the keypress)
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Confirm (Enter)', ...
          'Units', 'normalized', 'Position', [0.08 0.03 0.18 0.06], ...
          'FontSize', 10, 'FontWeight', 'bold', ...
          'BackgroundColor', [0.70 0.95 0.70], ...
          'Callback',     @(~,~) confirm_and_close(), ...
          'KeyPressFcn',  @on_key);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Reset shift (R)', ...
          'Units', 'normalized', 'Position', [0.29 0.03 0.16 0.06], ...
          'FontSize', 10, ...
          'Callback',     @(~,~) reset_and_redraw(), ...
          'KeyPressFcn',  @on_key);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Reset ROI (T)', ...
          'Units', 'normalized', 'Position', [0.48 0.03 0.16 0.06], ...
          'FontSize', 10, ...
          'Callback',     @(~,~) reset_roi(), ...
          'KeyPressFcn',  @on_key);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Skip (Esc)', ...
          'Units', 'normalized', 'Position', [0.67 0.03 0.18 0.06], ...
          'FontSize', 10, ...
          'BackgroundColor', [0.95 0.80 0.80], ...
          'Callback',     @(~,~) skip_and_close(), ...
          'KeyPressFcn',  @on_key);
update_overlay();
uiwait(fig);
% Save final diagnostic snapshot if possible
if isfield(params, 'log_dir') && ~isempty(params.log_dir) && (dr ~= 0 || dc ~= 0)
    fprintf('[ALIGN] Manual alignment committed: dr=%+d  dc=%+d px\n', dr, dc);
end
    % ------------- nested callbacks share dr/dc with the outer fn --------
    function update_overlay()
        if ~isgraphics(fig); return; end
        set(h_edge, 'XData', [1 + dc, W + dc], 'YData', [1 + dr, H + dr]);
        set(h_cent, 'XData', orig_cx + dc, 'YData', orig_cy + dr);
        for kk = 1:n_cent
            if isgraphics(h_txt(kk))
                set(h_txt(kk), 'Position', ...
                    [orig_cx(kk) + 8 + dc, orig_cy(kk) + dr, 0]);
            end
        end
        set(h_title, 'String', sprintf([ ...
            '10X + 40X GT overlay  |  Shift (dr=%+d, dc=%+d) px  =  (%+.1f, %+.1f) um\n' ...
            'Drag image (outside orange box) to shift GT  |  drag orange box to define scoring ROI\n' ...
            'Arrows: 1 px shift  |  Shift+arrow: 10 px  |  Enter: Confirm   Esc: Skip   R: Reset shift   T: Reset ROI'], ...
            dr, dc, dr * px_10x, dc * px_10x));
        drawnow limitrate;
    end
    function on_key(~, ev)
        step = 1;
        if ismember('shift', ev.Modifier); step = 10; end
        moved = false;
        switch ev.Key
            case 'uparrow';    dr = dr - step; moved = true;
            case 'downarrow';  dr = dr + step; moved = true;
            case 'leftarrow';  dc = dc - step; moved = true;
            case 'rightarrow'; dc = dc + step; moved = true;
            case 'return';     confirm_and_close(); return;
            case 'escape';     skip_and_close();    return;
            case 'r';          reset_and_redraw();  return;
            case 't';          reset_roi();         return;
        end
        if moved; update_overlay(); end
    end
    function on_button_down(~, ~)
        if ~isgraphics(ax); return; end
        pt = get(ax, 'CurrentPoint');
        xc = pt(1, 1); yc = pt(1, 2);
        % If the click lands inside (or near the handles of) the scoring
        % ROI, let the drawrectangle ROI handle it -- don't start an
        % alignment drag.  A small margin accommodates the corner/edge
        % drag handles that extend slightly outside h_roi.Position.
        if isgraphics(h_roi) && numel(h_roi.Position) == 4
            rp  = h_roi.Position;
            tol = 12;   % px of handle overhang
            if xc >= rp(1) - tol && xc <= rp(1) + rp(3) + tol && ...
               yc >= rp(2) - tol && yc <= rp(2) + rp(4) + tol
                return;
            end
        end
        if xc >= 0.5 && xc <= W + 0.5 && yc >= 0.5 && yc <= H + 0.5
            is_dragging   = true;
            drag_start_xy = [xc, yc];
            drag_start_dr = dr;
            drag_start_dc = dc;
            set(fig, 'Pointer', 'fleur');
        end
    end
    function on_button_motion(~, ~)
        if ~is_dragging || ~isgraphics(ax); return; end
        pt = get(ax, 'CurrentPoint');
        dx = pt(1, 1) - drag_start_xy(1);
        dy = pt(1, 2) - drag_start_xy(2);
        dc = round(drag_start_dc + dx);
        dr = round(drag_start_dr + dy);
        update_overlay();
    end
    function on_button_up(~, ~)
        if is_dragging
            is_dragging = false;
            set(fig, 'Pointer', 'arrow');
        end
    end
    function confirm_and_close()
        % Capture ROI position before the figure is torn down
        if isgraphics(h_roi) && numel(h_roi.Position) == 4
            pos = round(h_roi.Position);
            pos(1) = max(1, pos(1));
            pos(2) = max(1, pos(2));
            pos(3) = min(W - pos(1) + 1, max(1, pos(3)));
            pos(4) = min(H - pos(2) + 1, max(1, pos(4)));
            if pos(3) > 5 && pos(4) > 5
                roi_rect = pos;   % [x y w h]
            end
        end
        if isgraphics(fig); delete(fig); end
    end
    function skip_and_close()
        dr = 0; dc = 0; roi_rect = [];
        if isgraphics(fig); delete(fig); end
    end
    function reset_and_redraw()
        dr = 0; dc = 0;
        update_overlay();
    end
    function reset_roi()
        if isgraphics(h_roi)
            h_roi.Position = init_roi;
            drawnow;
        end
    end
end
function results = run_segmentation_validation(img_data, mask_path, csv_path, params)
%RUN_SEGMENTATION_VALIDATION  Compare segmentation methods vs annotated ground truth.
%
%  img_data   -- pre-loaded image struct from load_nd2 (cached in run_benchmark)
%  mask_path  -- multi-frame TIFF: Frame 0 = nucleus label mask (uint16)
%  csv_path   -- annotation CSV

fprintf('\n==============================================\n');
fprintf('  SEGMENTATION BENCHMARK vs GROUND TRUTH\n');
mip_str = '';
if isfield(params, 'gt_use_mip') && params.gt_use_mip; mip_str = ' | MIP'; end
fprintf('  FOV1  |  405nm DAPI channel  |  Series=%d%s\n', params.gt_series, mip_str);
fprintf('==============================================\n\n');

%% -- Step 1: Load ground truth mask and annotations -----------------------
fprintf('[BENCH] Loading ground truth...\n');

% Load nucleus mask (Frame 0)
gt_mask = load_gt_mask(mask_path, 0);
fprintf('[BENCH] GT mask: %dx%d, %d nuclei labelled\n', ...
        size(gt_mask,1), size(gt_mask,2), max(gt_mask(:)));

% Load annotations CSV (optional — ROI + morphology breakdown need it)
if ~isempty(csv_path) && isfile(csv_path)
    annot = readtable(csv_path);
    fprintf('[BENCH] Annotations: %d cells\n', height(annot));
    fprintf('[BENCH] Morphology breakdown:\n');
    morphs = unique(annot.morphology);
    for m = 1:numel(morphs)
        n_m = sum(strcmp(annot.morphology, morphs{m}));
        fprintf('          %-12s : %d\n', morphs{m}, n_m);
    end
else
    fprintf('[BENCH] No annotation CSV provided — ROI cropping disabled.\n');
    annot = table();   % empty; downstream ROI block will skip
    params.gt_roi_enabled = false;
end

%% -- Step 2: Use pre-loaded image data -----------------------------------
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

%% -- Step 2b: 10X/40X alignment + ROI selection (Option 7 pipeline) ------
% Runs the complete Option 7 alignment pipeline before any segmentation:
%   1. Build 10X stitched mosaic from stage coords
%   2. Load 40X FOVs, compute MIPs, read stage coords
%   3. Place 40X GT mask into 10X canvas (downsampled to 10X scale)
%   4. Crop overlap region (where both 10X and 40X exist)
%   5. Initial diagnostic plot
%   6. NCC refinement (shift GT to align with 10X segmentation)
%   7. Interactive ROI selection on the aligned 10X+GT overlap
%   8. Map confirmed 10X ROI -> 40X GT mask coords for downstream scoring
%
% If gt_roi_enabled=false the whole block is skipped.
% roi_rect = [x, y, w, h] in 40X GT mask pixel coords (applied below).

roi_rect = [];   % empty = score full image; [x y w h] in 40X gt_mask coords
roi_ids  = [];   % empty = no ID filter; otherwise list of GT cell IDs
                 %                        that fell inside the user's 10X ROI

if isfield(params, 'gt_roi_enabled') && params.gt_roi_enabled

    fprintf('\n[BENCH] Running 10X/40X alignment pipeline for ROI selection...\n');

    % -- 1. Build 10X mosaic -----------------------------------------------
    fprintf('[BENCH] Step 1: Building 10X mosaic...\n');
    [mosaic_10x, canvas_info, px_10x] = build_10x_mosaic_opt6(params);
    fprintf('[BENCH]   Mosaic: %dx%d px  (%.0fx%.0f um)\n', ...
            size(mosaic_10x,2), size(mosaic_10x,1), ...
            size(mosaic_10x,2)*px_10x, size(mosaic_10x,1)*px_10x);

    % -- 2. Load 40X FOVs + stage coords ----------------------------------
    fprintf('[BENCH] Step 2: Loading 40X FOVs and stage coords...\n');
    [~, fov_40x_coords, px_40x] = load_40x_mips_opt6(params);

    % -- 3. Place GT mask into 10X canvas ---------------------------------
    fprintf('[BENCH] Step 3: Placing GT mask into 10X canvas...\n');
    gt_canvas = place_gt_in_canvas_opt6( ...
        params, fov_40x_coords, canvas_info, px_10x, px_40x);

    % -- 4. Crop overlap region --------------------------------------------
    fprintf('[BENCH] Step 4: Cropping overlap region...\n');
    [crop_10x, crop_gt_c, ~] = crop_overlap_opt6(mosaic_10x, gt_canvas);
    fprintf('[BENCH]   Overlap: %dx%d px  (%.0fx%.0f um)  |  GT nuclei: %d\n', ...
            size(crop_10x,2), size(crop_10x,1), ...
            size(crop_10x,2)*px_10x, size(crop_10x,1)*px_10x, ...
            numel(unique(crop_gt_c(crop_gt_c>0))));

    % -- 5. Initial diagnostic ---------------------------------------------
    plot_overlap_diagnostic_opt6(crop_10x, crop_gt_c, px_10x, params, 'initial');

    % -- 6. NCC refinement -------------------------------------------------
    fprintf('[BENCH] Step 6: NCC alignment refinement...\n');
    seg_params_ncc             = params;
    seg_params_ncc.nucleus_channel = 1;
    seg_params_ncc.metadata.px_size_x = px_10x;
    [~, seg_initial] = segment_nuclei(crop_10x, seg_params_ncc);
    fprintf('[BENCH]   Initial 10X nuclei for NCC template: %d\n', ...
            numel(unique(seg_initial(seg_initial>0))));

    if ~isfield(params, 'align_ncc_search_um'); params.align_ncc_search_um = 100; end
    search_px = round(params.align_ncc_search_um / px_10x);
    [shift_r, shift_c, peak_ncc] = refine_alignment_ncc_opt6( ...
        seg_initial > 0, crop_gt_c > 0, search_px);
    fprintf('[BENCH]   Shift: (dr=%+d, dc=%+d) px = (%+.1f, %+.1f) um  NCC=%.4f\n', ...
            shift_r, shift_c, shift_r*px_10x, shift_c*px_10x, peak_ncc);

    if shift_r ~= 0 || shift_c ~= 0
        crop_gt_c = imtranslate(crop_gt_c, [shift_c, shift_r], ...
                                'FillValues', 0, 'Method', 'nearest');
        fprintf('[BENCH]   GT shifted. Post-shift GT nuclei: %d\n', ...
                numel(unique(crop_gt_c(crop_gt_c>0))));
        plot_overlap_diagnostic_opt6(crop_10x, crop_gt_c, px_10x, params, 'refined');
    end

    % -- 7. Interactive ROI selection (identical to Option 7 Step 5.7) ----
    roi_rect_10x = [];
    roi_enabled  = true;

    if roi_enabled
        fprintf('\n[BENCH] Step 7: ROI selection on aligned 10X+GT overlap...\n');
        fprintf('[BENCH]   Draw a rectangle, then DOUBLE-CLICK inside it.\n');
        fprintf('[BENCH]   Close the window to skip and score the full crop.\n\n');

        gt_bw_sel   = crop_gt_c > 0;
        gt_edge_sel = edge(gt_bw_sel, 'canny');
        roi_rgb     = repmat(mat2gray(crop_10x), [1 1 3]);
        roi_rgb(:,:,1) = min(1, roi_rgb(:,:,1) + gt_edge_sel);

        gt_ids_sel   = unique(crop_gt_c(crop_gt_c > 0));
        gt_cents_sel = regionprops(crop_gt_c, 'Centroid');
        [Hcrop, Wcrop] = size(crop_10x);

        roi_fig = figure('Name', ...
            'Option 6 -- ROI Selection on 10X FOV (draw on 10X image, double-click to confirm)', ...
            'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.04 0.04 0.72 0.88]);
        ax_sel = axes(roi_fig, 'Units', 'normalized', 'Position', [0.05 0.08 0.90 0.82]);
        imagesc(ax_sel, roi_rgb);
        axis(ax_sel, 'image'); axis(ax_sel, 'off');
        hold(ax_sel, 'on');

        for k = 1:numel(gt_ids_sel)
            id = gt_ids_sel(k);
            if id <= numel(gt_cents_sel) && ~isempty(gt_cents_sel(id).Centroid)
                cx = gt_cents_sel(id).Centroid(1);
                cy = gt_cents_sel(id).Centroid(2);
                plot(ax_sel, cx, cy, 'co', 'MarkerSize', 9, 'LineWidth', 1.5, ...
                     'MarkerFaceColor', 'none');
                text(ax_sel, cx + 8, cy, num2str(id), ...
                     'Color', 'cyan', 'FontSize', 6, 'FontWeight', 'bold');
            end
        end

        plot(ax_sel, [1 Wcrop Wcrop 1 1], [1 1 Hcrop Hcrop 1], 'y--', 'LineWidth', 1.5);
        text(ax_sel, Wcrop*0.02, Hcrop*0.04, ...
             sprintf('Full crop: %d x %d px  (%.0f x %.0f um)', ...
                     Wcrop, Hcrop, Wcrop*px_10x, Hcrop*px_10x), ...
             'Color', 'yellow', 'FontSize', 8, 'BackgroundColor', [0 0 0 0.45]);

        title(ax_sel, ...
            {['10X Mosaic Overlap  |  NCC-refined  |  ' ...
              sprintf('%d GT nuclei (cyan circles + red boundaries)', numel(gt_ids_sel))], ...
             'Draw a rectangle to define the valid scoring zone, then DOUBLE-CLICK inside it', ...
             'Close this window to skip ROI and score the full crop'}, ...
            'FontSize', 10, 'FontWeight', 'bold');

        try
            h_rect = drawrectangle(ax_sel, ...
                         'Color',          [1.0 0.55 0.0], ...
                         'LineWidth',       2.5, ...
                         'Label',          'Valid Scoring ROI -- double-click to confirm', ...
                         'LabelTextColor', [1 1 1], ...
                         'FaceAlpha',       0.07);
            wait(h_rect);
            roi_rect_10x = round(h_rect.Position);
            roi_rect_10x(1) = max(1, roi_rect_10x(1));
            roi_rect_10x(2) = max(1, roi_rect_10x(2));
            roi_rect_10x(3) = min(roi_rect_10x(3), Wcrop - roi_rect_10x(1) + 1);
            roi_rect_10x(4) = min(roi_rect_10x(4), Hcrop - roi_rect_10x(2) + 1);

            fprintf('[BENCH] ROI (10X): x=%d  y=%d  w=%d  h=%d  (%.0f x %.0f um)\n', ...
                    roi_rect_10x(1), roi_rect_10x(2), ...
                    roi_rect_10x(3), roi_rect_10x(4), ...
                    roi_rect_10x(3)*px_10x, roi_rect_10x(4)*px_10x);

            rectangle(ax_sel, 'Position', roi_rect_10x, ...
                      'EdgeColor', [1 0.55 0], 'LineWidth', 3.0);
            title(ax_sel, ...
                {sprintf('ROI confirmed: x=%d  y=%d  w=%d  h=%d  (%.0f x %.0f um)', ...
                         roi_rect_10x(1), roi_rect_10x(2), ...
                         roi_rect_10x(3), roi_rect_10x(4), ...
                         roi_rect_10x(3)*px_10x, roi_rect_10x(4)*px_10x), ...
                 'Mapping to 40X GT space...'}, ...
                'FontSize', 10, 'FontWeight', 'bold');
            drawnow;

        catch ME_roi
            warning('[BENCH] ROI draw cancelled (%s). Scoring full image.', ME_roi.message);
            roi_rect_10x = [];
        end

        if isfield(params, 'log_dir') && ~isempty(params.log_dir)
            exportgraphics(roi_fig, fullfile(params.log_dir, ...
                sprintf('roi_selection_%s.png', params.run_timestamp)), 'Resolution', 150);
            fprintf('[BENCH] ROI figure saved.\n');
        end
        pause(1.5);
        close(roi_fig);
    end

    % -- 8. Map 10X ROI -> 40X: label lookup + direct rectangle mapping ----
    % Step (a): identify the GT cell IDs that appear under the user's
    %           rectangle on the 10X overlay.  imresize with 'nearest'
    %           preserves label values, so 10X labels == 40X labels.
    % Step (b): map the drawn rectangle itself from 10X overlay coords to
    %           40X gt_mask coords so the physical crop matches exactly
    %           what the user outlined.  We derive the affine transform
    %           empirically by using one cell in roi_ids as a reference:
    %           its centroid is known in both frames, so we can solve for
    %           the offset around a uniform scale = px_10x / px_40x.
    if ~isempty(roi_rect_10x)
        [Hcc, Wcc] = size(crop_gt_c);
        rx_c = max(1, roi_rect_10x(1));
        ry_c = max(1, roi_rect_10x(2));
        rw_c = min(Wcc - rx_c + 1, roi_rect_10x(3));
        rh_c = min(Hcc - ry_c + 1, roi_rect_10x(4));

        % (a) label lookup
        roi_patch = crop_gt_c(ry_c:(ry_c+rh_c-1), rx_c:(rx_c+rw_c-1));
        roi_ids   = unique(roi_patch(roi_patch > 0));
        fprintf('[BENCH] User ROI covers %d GT cell IDs\n', numel(roi_ids));

        if ~isempty(roi_ids)
            % (b) derive transform with one reference cell in roi_ids
            scale_10_to_40 = px_10x / px_40x;
            offset_xy      = [];
            for kref = 1:numel(roi_ids)
                ref_id   = roi_ids(kref);
                stats_cc = regionprops(crop_gt_c == ref_id, 'Centroid');
                stats_gm = regionprops(gt_mask   == ref_id, 'Centroid');
                if ~isempty(stats_cc) && ~isempty(stats_gm) ...
                        && ~isempty(stats_cc(1).Centroid) ...
                        && ~isempty(stats_gm(1).Centroid)
                    cent_cc    = stats_cc(1).Centroid;    % [x, y] in 10X overlay
                    cent_gm    = stats_gm(1).Centroid;    % [x, y] in 40X gt_mask
                    offset_xy  = cent_gm - cent_cc * scale_10_to_40;
                    break;
                end
            end

            if ~isempty(offset_xy)
                roi_40x_x = round(rx_c   * scale_10_to_40 + offset_xy(1));
                roi_40x_y = round(ry_c   * scale_10_to_40 + offset_xy(2));
                roi_40x_w = round(rw_c   * scale_10_to_40);
                roi_40x_h = round(rh_c   * scale_10_to_40);

                [Hgt, Wgt] = size(gt_mask);
                roi_40x_x  = max(1, roi_40x_x);
                roi_40x_y  = max(1, roi_40x_y);
                roi_40x_w  = min(roi_40x_w, Wgt - roi_40x_x + 1);
                roi_40x_h  = min(roi_40x_h, Hgt - roi_40x_y + 1);

                roi_rect = [roi_40x_x, roi_40x_y, roi_40x_w, roi_40x_h];
                fprintf('[BENCH] ROI mapped 10X->40X (rect, not bbox): x=%d y=%d w=%d h=%d\n', ...
                        roi_rect(1), roi_rect(2), roi_rect(3), roi_rect(4));
                fprintf('[BENCH]   transform: scale=%.4f  offset=(%.1f, %.1f)  ref_cell=%d\n', ...
                        scale_10_to_40, offset_xy(1), offset_xy(2), ref_id);
            else
                % Fallback: tight bbox of cells (no valid reference centroid found)
                in_roi_pix = ismember(gt_mask, roi_ids);
                [rr, cc]   = find(in_roi_pix);
                if ~isempty(rr)
                    x0 = min(cc);  y0 = min(rr);
                    roi_rect = [x0, y0, max(cc) - x0 + 1, max(rr) - y0 + 1];
                    fprintf('[BENCH] No valid reference cell -- falling back to cells bbox.\n');
                end
            end
        else
            fprintf('[BENCH] ROI contains no GT cells -- scoring full image.\n');
        end
    end

else
    fprintf('[BENCH] ROI selection disabled -- scoring full image.\n');
end

% -- Apply ROI crop to image, GT mask, and annotation table ---------------
% The scoring set after this block is exactly the GT cells whose labels
% appeared in the user's drawn rectangle on the 10X overlay (roi_ids).
% We use label-based filtering (ismember on cell IDs) rather than a
% centroid-in-rectangle test, because the drawn rectangle is defined in
% the 10X overlay's coordinate frame and cells can straddle its edges.
% Original GT cell IDs are preserved -- no renumbering -- so
% plot_benchmark_results' `find(gt_ids == annot.cell_id)` lookup works.
if ~isempty(roi_rect)
    rx = roi_rect(1); ry = roi_rect(2);
    rw = roi_rect(3); rh = roi_rect(4);
    row_range = ry : (ry + rh - 1);
    col_range = rx : (rx + rw - 1);

    fprintf('[BENCH] Applying ROI crop...\n');

    % 1) Filter annotation.  Apply BOTH: label membership in roi_ids (the
    %    semantic selection from the 10X overlay) AND centroid-in-rect
    %    (so cells whose centroid would fall outside the physical 40X
    %    crop don't become phantom entries after the crop zeros their
    %    labels out).  The intersection guarantees annot and the cropped
    %    gt_mask share the same cell IDs.
    n_annot_before = height(annot);
    in_rect = annot.centroid_x >= rx & annot.centroid_x <= (rx + rw - 1) & ...
              annot.centroid_y >= ry & annot.centroid_y <= (ry + rh - 1);
    if ~isempty(roi_ids)
        annot = annot(ismember(annot.cell_id, roi_ids) & in_rect, :);
    else
        annot = annot(in_rect, :);
    end
    fprintf('[BENCH]   Annotations after ROI filter: %d / %d cells retained\n', ...
            height(annot), n_annot_before);

    % 2) Zero out any GT mask pixel whose label is NOT in the final annot
    %    cell list.  That way the cropped gt_mask contains exactly the
    %    cells displayed by plot_benchmark_results -- no orphan labels.
    n_gt_before = numel(unique(gt_mask(gt_mask > 0)));
    if ~isempty(annot)
        keep_pixels = ismember(gt_mask, annot.cell_id);
    else
        keep_pixels = false(size(gt_mask));
    end
    gt_mask(~keep_pixels) = 0;
    n_gt_after = numel(unique(gt_mask(gt_mask > 0)));
    fprintf('[BENCH]   GT mask labels after ROI filter: %d / %d kept (original IDs preserved)\n', ...
            n_gt_after, n_gt_before);

    % 3) Physically crop image + GT mask to ROI bounding box.
    nucleus_img = nucleus_img(row_range, col_range);
    gt_mask     = gt_mask(row_range, col_range);

    % 4) Remap annotation centroid coords to the cropped-image frame.
    annot.centroid_x = annot.centroid_x - (rx - 1);
    annot.centroid_y = annot.centroid_y - (ry - 1);

    % 5) Update img_data metadata to reflect the cropped dimensions.
    img_data.metadata.img_width  = rw;
    img_data.metadata.img_height = rh;
end

%% -- Step 3: Run each segmentation method --------------------------------
fprintf('\n[BENCH] Running segmentation methods...\n');

% Method 1: Otsu baseline (L1)
fprintf('[BENCH] Method 1: Otsu baseline (L1)...\n');
[~, seg_matlab] = segment_nuclei_l1(nucleus_img, params);
fprintf('[BENCH]   Detected %d connected components\n', max(seg_matlab(:)));

% Method 2: Rule-based pipeline + circularity filter (L2)
% Uses segment_nuclei (adaptive threshold + morphology + watershed) as its
% base, then applies the circularity rule.
fprintf('[BENCH] Method 2: Rule-based + circularity filter (L2)...\n');
[cells_rule, seg_rule] = segment_nuclei(nucleus_img, params);
circ_thresh_val = min(params.circularity_threshold, 0.60);
circ_vals    = [cells_rule.circularity];
pass_mask_L2 = circ_vals >= circ_thresh_val;
seg_L2       = zeros(size(seg_rule), 'uint16');
cell_ids_L2  = find(pass_mask_L2);
for k = 1:numel(cell_ids_L2)
    seg_L2(seg_rule == cell_ids_L2(k)) = k;
end
fprintf('[BENCH]   thresh=%.2f  %d/%d cells pass circularity\n', ...
        circ_thresh_val, sum(pass_mask_L2), numel(cells_rule));

% Method 3: StarDist via bridge
fprintf('[BENCH] Method 3: StarDist (L3)...\n');
seg_stardist = run_stardist_full_image(nucleus_img, params);
n_sd = max(seg_stardist(:));
fprintf('[BENCH]   StarDist detected %d nuclei\n', n_sd);

%% -- Step 4: Compute IoU-based ROC for each method -----------------------
fprintf('\n[BENCH] Computing ROC curves...\n');

methods = {'Otsu (L1)', 'Rule-based + Circ. (L2)', 'StarDist (L3)'};
segs    = {seg_matlab, seg_L2, seg_stardist};
colors  = {[0.3 0.6 1.0], [1.0 0.4 0.3], [0.2 0.8 0.3]};
roc     = cell(1, 3);

for m = 1:3
    roc{m} = compute_roc(gt_mask, segs{m}, annot);
    P = roc{m}.tp_at_half / max(1, roc{m}.n_pred);
    R = roc{m}.tp_at_half / max(1, roc{m}.n_gt);
    roc{m}.f1_iou05 = 2*P*R / (P + R + eps);
    fprintf('[BENCH] %s -- AP=%.3f  F1@IoU=0.5=%.3f  TP=%d/%d (rec=%.2f)  pred=%d (prec=%.2f)\n', ...
            methods{m}, roc{m}.ap, roc{m}.f1_iou05, ...
            roc{m}.tp_at_half, roc{m}.n_gt, R, roc{m}.n_pred, P);
end

%% -- Step 5: Display results ----------------------------------------------
plot_benchmark_results(nucleus_img, gt_mask, segs, roc, methods, colors, annot, params, roi_rect);

%% -- Step 6 (optional): Export centroid CSV + accuracy histogram ----------
export_on = isfield(params, 'gt_export_centroids_csv') && params.gt_export_centroids_csv;
if export_on
    fprintf('\n[BENCH] Exporting centroid CSVs and accuracy histogram...\n');

    % Pixel size for converting px -> um coordinates
    px_um_x = img_data.metadata.px_size_x;
    px_um_y = img_data.metadata.px_size_y;

    method_tags = {'L1', 'L2', 'L3'};
    acc_data    = struct();   % accumulate per-method accuracy for histogram

    for m = 1:3
        seg       = segs{m};
        pred_ids  = unique(seg(seg > 0));
        n_pred    = numel(pred_ids);
        iou_mat   = roc{m}.iou_matrix;
        gt_ids    = unique(gt_mask(gt_mask > 0));
        n_gt      = numel(gt_ids);

        % -- Greedy IoU matching at threshold 0.5 -------------------------
        matched_gt   = false(n_gt, 1);
        matched_pred = false(n_pred, 1);
        if n_pred > 0 && size(iou_mat,1) == n_gt && size(iou_mat,2) == n_pred
            [sorted_iou, sort_idx] = sort(iou_mat(:), 'descend');
            for si = 1:numel(sorted_iou)
                if sorted_iou(si) < 0.5; break; end
                [gi, pi] = ind2sub([n_gt, n_pred], sort_idx(si));
                if ~matched_gt(gi) && ~matched_pred(pi)
                    matched_gt(gi)   = true;
                    matched_pred(pi) = true;
                end
            end
        end

        % -- Build centroid table for this method -------------------------
        cell_id_out  = zeros(n_pred, 1);
        centroid_x_px= zeros(n_pred, 1);
        centroid_y_px= zeros(n_pred, 1);
        centroid_x_um= zeros(n_pred, 1);
        centroid_y_um= zeros(n_pred, 1);
        matched_gt_id= zeros(n_pred, 1);   % 0 = no GT match
        iou_score    = zeros(n_pred, 1);

        for pi = 1:n_pred
            mask_p = (seg == pred_ids(pi));
            props  = regionprops(mask_p, 'Centroid');
            if ~isempty(props)
                cx = props(1).Centroid(1);
                cy = props(1).Centroid(2);
            else
                cx = 0; cy = 0;
            end
            cell_id_out(pi)   = pred_ids(pi);
            centroid_x_px(pi) = cx;
            centroid_y_px(pi) = cy;
            centroid_x_um(pi) = cx * px_um_x;
            centroid_y_um(pi) = cy * px_um_y;

            % Best GT match for this prediction
            if size(iou_mat,2) >= pi
                [best_iou, best_gi] = max(iou_mat(:, pi));
                iou_score(pi)       = best_iou;
                if matched_pred(pi)
                    matched_gt_id(pi) = gt_ids(best_gi);
                end
            end
        end

        % -- Write CSV -----------------------------------------------------
        T = table(cell_id_out, centroid_x_px, centroid_y_px, ...
                  centroid_x_um, centroid_y_um, ...
                  matched_gt_id, iou_score, ...
                  'VariableNames', {'cell_id', ...
                                    'centroid_x_px', 'centroid_y_px', ...
                                    'centroid_x_um', 'centroid_y_um', ...
                                    'matched_gt_id', 'iou_score'});

        csv_out = fullfile(params.log_dir, ...
            sprintf('centroids_%s_%s.csv', method_tags{m}, params.run_timestamp));
        writetable(T, csv_out);
        fprintf('[BENCH] %s: centroid CSV saved (%d cells) -> %s\n', ...
                method_tags{m}, n_pred, csv_out);

        % -- Accuracy buckets (% of GT detected at increasing IoU thresholds)
        iou_bins = 0:0.05:1.0;   % x-axis buckets: IoU >= bin value
        det_pct  = zeros(size(iou_bins));
        for bi = 1:numel(iou_bins)
            % Count GT cells that have at least one prediction with IoU >= bin
            if n_pred > 0 && size(iou_mat,2) == n_pred
                best_iou_per_gt = max(iou_mat, [], 2);  % [n_gt x 1]
                det_pct(bi) = sum(best_iou_per_gt >= iou_bins(bi)) / n_gt * 100;
            end
        end
        acc_data(m).label   = methods{m};
        acc_data(m).tag     = method_tags{m};
        acc_data(m).iou_bins= iou_bins;
        acc_data(m).det_pct = det_pct;
        acc_data(m).n_gt    = n_gt;
        acc_data(m).n_pred  = n_pred;
        acc_data(m).tp      = sum(matched_pred);
        acc_data(m).fp      = sum(~matched_pred);
        acc_data(m).fn      = sum(~matched_gt);
    end

    % -- Plot accuracy histogram figure ------------------------------------
    plot_accuracy_histogram(acc_data, roc, methods, colors, params);
end

%% -- Package output -------------------------------------------------------
results.roc        = roc;
results.methods    = methods;
results.gt_mask    = gt_mask;
results.segs       = segs;
results.annot      = annot;
results.nucleus_img= nucleus_img;
results.roi_rect   = roi_rect;   % [] = full image; [x y w h] = cropped ROI
if export_on
    results.acc_data = acc_data;
end
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
%COMPUTE_ROC  Compute ROC curve by varying confidence score threshold.
%
%  Confidence score = mean intensity of each predicted nucleus region.
%  At each score threshold:
%    - Keep only predictions with score >= threshold
%    - Match kept predictions to GT nuclei at fixed IoU=0.5
%    - Compute TPR = TP/n_GT and FPR = FP/(FP+n_GT)
%  This sweeps the full [0,1] range like a standard detector ROC.

iou_fixed = 0.5;   % fixed IoU threshold for matching

gt_ids   = unique(gt_mask(gt_mask > 0));
pred_ids = unique(pred_mask(pred_mask > 0));
n_gt     = numel(gt_ids);
n_pred   = numel(pred_ids);

fprintf('    Computing IoU matrix (%d GT x %d pred)...\n', n_gt, n_pred);

% Compute IoU matrix at fixed threshold
iou_matrix = zeros(n_gt, n_pred);
pred_scores = zeros(1, n_pred);

for pi = 1:n_pred
    pred_mask_p = (pred_mask == pred_ids(pi));
    % Score = mean intensity proxy (use mask area as fallback)
    pred_scores(pi) = sum(pred_mask_p(:));
    for gi = 1:n_gt
        gt_mask_i    = (gt_mask == gt_ids(gi));
        intersection = sum(gt_mask_i(:) & pred_mask_p(:));
        if intersection > 0
            union = sum(gt_mask_i(:) | pred_mask_p(:));
            iou_matrix(gi, pi) = intersection / union;
        end
    end
end

% Normalise scores to [0,1]
if max(pred_scores) > 0
    pred_scores = pred_scores ./ max(pred_scores);
end

% Sweep score thresholds from high to low (like a detector ROC)
% At threshold=1: keep nothing -> TPR=0, FPR=0
% At threshold=0: keep everything -> TPR=max, FPR=max
score_thresholds = [1.0, sort(unique(pred_scores), 'descend'), -eps];
n_thresh = numel(score_thresholds);
tpr      = zeros(1, n_thresh);
fpr      = zeros(1, n_thresh);
precision= zeros(1, n_thresh);

for ti = 1:n_thresh
    st = score_thresholds(ti);
    keep = pred_scores >= st;   % predictions above score threshold
    keep_idx = find(keep);
    n_keep   = numel(keep_idx);

    matched_gt   = false(n_gt, 1);
    matched_pred = false(n_keep, 1);

    if n_keep > 0
        % Greedy match at IoU=0.5
        sub_iou = iou_matrix(:, keep_idx);
        [sorted_iou, sort_idx] = sort(sub_iou(:), 'descend');
        for si = 1:numel(sorted_iou)
            if sorted_iou(si) < iou_fixed; break; end
            [gi, pi] = ind2sub([n_gt, n_keep], sort_idx(si));
            if ~matched_gt(gi) && ~matched_pred(pi)
                matched_gt(gi)   = true;
                matched_pred(pi) = true;
            end
        end
    end

    TP = sum(matched_gt);
    FP = n_keep - sum(matched_pred);

    tpr(ti)  = TP / (n_gt + eps);
    fpr(ti)  = FP / (FP + n_gt + eps);
    if TP + FP > 0
        precision(ti) = TP / (TP + FP);
    end
end

% AUC via trapezoidal integration (sort by FPR)
[fpr_s, idx] = sort(fpr);
tpr_s        = tpr(idx);
auc          = trapz(fpr_s, tpr_s);

% F1 at IoU=0.5, best threshold
t50_idx    = find(thresholds_f1(tpr, precision), 1);
if isempty(t50_idx); t50_idx = 1; end
% Use score threshold closest to 0.5 for F1 reference
[~, mid]   = min(abs(score_thresholds - 0.5));
f1_at_half = 2 * precision(mid) * tpr(mid) / (precision(mid) + tpr(mid) + eps);

roc.thresholds  = score_thresholds;
roc.tpr         = tpr;
roc.fpr         = fpr;
roc.precision   = precision;
roc.auc         = auc;
roc.f1_at_half  = f1_at_half;
roc.n_gt        = n_gt;
roc.n_pred      = n_pred;
roc.iou_matrix  = iou_matrix;
roc.pred_scores = pred_scores;

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
    [si_v, idx_v] = sort(iou_matrix(:), 'descend');
    for k = 1:numel(si_v)
        if si_v(k) < thr; break; end
        [gi, pj] = ind2sub([n_gt, n_pred], idx_v(k));
        if ~m_gt(gi) && ~m_pred(pj)
            m_gt(gi) = true; m_pred(pj) = true;
        end
    end
    TP = sum(m_gt);
    if n_pred > 0; pr_prec(ii) = TP / n_pred; else; pr_prec(ii) = 1; end
    if n_gt   > 0; pr_rec(ii)  = TP / n_gt;   else; pr_rec(ii)  = 0; end
    if abs(thr - 0.5) < 1e-9; tp_at_half = TP; end
end
[rec_s, sord] = sort(pr_rec);
prec_s        = pr_prec(sord);
ap            = trapz([0 rec_s 1], [prec_s(1) prec_s prec_s(end)]);

roc.iou_grid     = iou_grid;
roc.pr_precision = pr_prec;
roc.pr_recall    = pr_rec;
roc.ap           = ap;
roc.tp_at_half   = tp_at_half;
end

function idx = thresholds_f1(tpr, precision)
% Helper: find index of best F1
f1 = 2 * precision .* tpr ./ (precision + tpr + eps);
[~, idx] = max(f1);
end


%% ======================================================================
function plot_benchmark_results(nucleus_img, gt_mask, segs, roc, methods, colors, annot, params, roi_rect)
%PLOT_BENCHMARK_RESULTS  Option-7-style 2x3 layout: method overlays,
%  ROC curves, detection counts, and a 10X+GT reference panel.
%
%  Method panels show predicted boundaries in the method's color plus GT
%  boundaries in white -- no +/X markers, no morphology circles.  The
%  subtitle reports AUC, F1, and number of predicted objects.

if nargin < 9; roi_rect = []; end

roi_str = ' (full crop)';
if ~isempty(roi_rect)
    roi_str = sprintf(' (ROI: x=%d y=%d w=%d h=%d)', ...
                      roi_rect(1), roi_rect(2), roi_rect(3), roi_rect(4));
end

fig = figure('Name', 'Segmentation Benchmark vs Ground Truth', ...
             'NumberTitle', 'off', 'Units', 'normalized', ...
             'Position', [0.02 0.02 0.96 0.92]);
sgtitle(['Segmentation Validation -- 40X input vs 40X GT mask' roi_str], ...
        'FontSize', 12, 'FontWeight', 'bold');

method_short = {'Otsu (L1)', 'Rule-based (L2)', 'StarDist (L3)'};

%% Panels 1-3: method overlays (predicted edges + GT edges)
for m = 1:3
    ax = subplot(2, 3, m);
    pred_edge = edge(segs{m} > 0, 'canny');
    gt_edge   = edge(gt_mask > 0, 'canny');
    rgb = repmat(mat2gray(nucleus_img), [1 1 3]);
    rgb(:,:,1) = rgb(:,:,1) + pred_edge * colors{m}(1) * 0.8;
    rgb(:,:,2) = rgb(:,:,2) + pred_edge * colors{m}(2) * 0.8;
    rgb(:,:,3) = rgb(:,:,3) + pred_edge * colors{m}(3) * 0.8;
    rgb(:,:,:) = rgb(:,:,:) + repmat(gt_edge * 0.5, [1 1 3]);
    imagesc(ax, min(1, rgb)); axis(ax,'image'); axis(ax,'off');

    P_m = roc{m}.tp_at_half / max(1, roc{m}.n_pred);
    R_m = roc{m}.tp_at_half / max(1, roc{m}.n_gt);
    f1_iou05 = 2*P_m*R_m / (P_m + R_m + eps);
    title(ax, sprintf('%s\nAP=%.3f  F1@IoU=0.5=%.3f  matched=%d/%d', ...
          method_short{m}, roc{m}.ap, f1_iou05, ...
          roc{m}.tp_at_half, roc{m}.n_gt), ...
          'FontSize', 9);
end

%% Panel 4: Precision-Recall curves over IoU threshold sweep (COCO-style)
ax4 = subplot(2, 3, 4);
hold(ax4, 'on');
for m = 1:3
    [rec_s, si] = sort(roc{m}.pr_recall);
    prec_s      = roc{m}.pr_precision(si);
    plot(ax4, rec_s, prec_s, '-o', 'Color', colors{m}, 'LineWidth', 2, ...
         'MarkerSize', 4, 'MarkerFaceColor', colors{m}, ...
         'DisplayName', sprintf('%s (AP=%.3f)', method_short{m}, roc{m}.ap));
    [~, i05] = min(abs(roc{m}.iou_grid - 0.5));
    plot(ax4, roc{m}.pr_recall(i05), roc{m}.pr_precision(i05), ...
         'p', 'MarkerFaceColor', colors{m}, 'MarkerEdgeColor', 'k', ...
         'MarkerSize', 10, 'HandleVisibility', 'off');
end
xlabel(ax4, 'Recall');
ylabel(ax4, 'Precision');
title(ax4, 'Precision-Recall (IoU sweep 0.05-0.95)', 'FontSize', 10);
legend(ax4, 'Location', 'southwest', 'FontSize', 8);
xlim(ax4, [0 1]); ylim(ax4, [0 1.05]);
grid(ax4, 'on'); axis(ax4, 'square');
text(ax4, 0.02, 0.08, 'star = IoU 0.5', 'FontSize', 7, 'Color', [0.3 0.3 0.3]);

%% Panel 5: Detection Rate (% of GT cells)
ax5 = subplot(2, 3, 5);
bar_data = zeros(3, 3);
for m = 1:3
    n_gt_m = max(1, roc{m}.n_gt);
    bar_data(m, 1) = 100;                                  % GT = 100% baseline
    bar_data(m, 2) = 100 * roc{m}.n_pred    / n_gt_m;      % predicted as % of GT
    bar_data(m, 3) = 100 * roc{m}.tp_at_half / n_gt_m;     % matched @ IoU=0.5 (recall %)
end
b = bar(ax5, bar_data');
for m = 1:3; b(m).FaceColor = colors{m}; end
set(ax5, 'XTickLabel', {'GT nuclei','Predicted','Matched @ IoU=0.5'}, 'XTick', 1:3);
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
grid(ax5, 'on');

%% Panel 6: Input + GT boundaries reference
ax6 = subplot(2, 3, 6);
rgb = repmat(mat2gray(nucleus_img), [1 1 3]);
gt_edge = edge(gt_mask > 0, 'canny');
rgb(:,:,1) = min(1, rgb(:,:,1) + gt_edge);
imagesc(ax6, rgb); axis(ax6, 'image'); axis(ax6, 'off');
title(ax6, '40X + GT boundaries', 'FontSize', 10);

drawnow;

% Save
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_path = fullfile(params.log_dir, ...
        sprintf('segmentation_benchmark_%s.png', params.run_timestamp));
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Benchmark figure saved.\n');
end
end


%% ======================================================================
function plot_accuracy_histogram(acc_data, roc, methods, colors, params)
%PLOT_ACCURACY_HISTOGRAM  Three paired bar panels: L1/GT, L2/GT, L3/GT
%
%  Each panel shows detection rate (%) vs IoU threshold sweep.
%  The GT reference bar (100%) is drawn as a horizontal dashed line.
%  A summary panel (Panel 4) overlays all three methods for comparison.
%
%  Layout:  [L1/GT] [L2/GT] [L3/GT] [All-method overlay]

fig = figure('Name', 'Segmentation Accuracy vs Ground Truth', ...
             'NumberTitle', 'off', 'Units', 'normalized', ...
             'Position', [0.04 0.08 0.92 0.82]);
sgtitle({'Segmentation Accuracy vs Ground Truth', ...
         'Y-axis: % of GT cells detected at IoU \geq threshold  |  X-axis: IoU threshold'}, ...
        'FontSize', 13, 'FontWeight', 'bold');

method_short = {'Otsu (L1)', 'Rule-based (L2)', 'StarDist (L3)'};
pair_titles  = {'L1 / GT', 'L2 / GT', 'L3 / GT'};

% Bar display at discrete IoU steps for clean histogram look
iou_display = 0:0.1:1.0;   % 11 bars per panel

% -- Panels 1-3: one per method vs GT -------------------------------------
for m = 1:3
    ax = subplot(2, 4, m);

    % Interpolate acc_data to display bins
    det_pct_disp = interp1(acc_data(m).iou_bins, acc_data(m).det_pct, ...
                           iou_display, 'linear', 0);
    det_pct_disp = max(0, det_pct_disp);

    % GT reference: 100% at IoU = 0 (all cells detectable at zero threshold)
    % Draw method bars
    bar(ax, iou_display, det_pct_disp, 0.7, ...
        'FaceColor', colors{m}, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    hold(ax, 'on');

    % GT = 100% reference line
    yline(ax, 100, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.8, ...
          'Label', 'GT (100%)', 'LabelHorizontalAlignment', 'left', ...
          'FontSize', 8);

    % Mark IoU = 0.5 operating point
    det_at_half = interp1(acc_data(m).iou_bins, acc_data(m).det_pct, 0.5, 'linear', 0);
    plot(ax, 0.5, det_at_half, 'kv', 'MarkerSize', 9, 'MarkerFaceColor', colors{m}, ...
         'LineWidth', 1.5);
    text(ax, 0.52, det_at_half + 3, sprintf('%.1f%%', det_at_half), ...
         'FontSize', 8, 'FontWeight', 'bold', 'Color', colors{m} * 0.7);

    % Precision / Recall annotation box
    tp  = acc_data(m).tp;
    fp  = acc_data(m).fp;
    fn  = acc_data(m).fn;
    prec = tp / (tp + fp + eps) * 100;
    rec  = tp / (tp + fn + eps) * 100;
    f1   = 2 * prec * rec / (prec + rec + eps);
    ann_str = sprintf('TP=%d  FP=%d  FN=%d\nPrec=%.1f%%  Rec=%.1f%%\nF1=%.1f%%  AUC=%.3f', ...
                      tp, fp, fn, prec, rec, f1, roc{m}.auc);
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

% -- Panel 4 (top-right): overlay of all three methods --------------------
ax4 = subplot(2, 4, 4);
hold(ax4, 'on');
leg_h = gobjects(3, 1);
for m = 1:3
    det_pct_full = acc_data(m).det_pct;
    iou_full     = acc_data(m).iou_bins;
    leg_h(m) = plot(ax4, iou_full, det_pct_full, '-o', ...
                    'Color', colors{m}, 'LineWidth', 2.2, 'MarkerSize', 5, ...
                    'DisplayName', method_short{m});
end
yline(ax4, 100, 'k--', 'LineWidth', 1.5, 'Label', 'GT (100%)', ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 8);
xline(ax4, 0.5, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, ...
      'Label', 'IoU=0.5', 'LabelVerticalAlignment', 'bottom', 'FontSize', 8);
xlabel(ax4, 'IoU Threshold');
ylabel(ax4, 'Detection Rate (% of GT)');
legend(ax4, leg_h, 'Location', 'southwest', 'FontSize', 8);
title(ax4, {'All Methods -- Overlay', 'detection rate vs IoU threshold'}, ...
      'FontSize', 10, 'FontWeight', 'bold');
xlim(ax4, [0 1]); ylim(ax4, [0 115]);
grid(ax4, 'on');

% -- Panel 5-6 (bottom left, spanning 2 cols): IoU=0.5 paired bar chart ---
ax5 = subplot(2, 4, [5 6]);
bar_vals = zeros(3, 1);
for m = 1:3
    bar_vals(m) = interp1(acc_data(m).iou_bins, acc_data(m).det_pct, 0.5, 'linear', 0);
end
b = bar(ax5, 1:3, bar_vals, 0.55, 'EdgeColor', 'none');
for m = 1:3
    b.FaceColor = 'flat';
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
set(ax5, 'XTick', 1:3, 'XTickLabel', {'L1 / GT', 'L2 / GT', 'L3 / GT'}, ...
         'FontSize', 10);
ylabel(ax5, 'Detection Rate at IoU \geq 0.5 (%)');
title(ax5, 'Paired Accuracy at IoU = 0.5  (L1/GT  +/-  L2/GT  +/-  L3/GT)', ...
      'FontSize', 11, 'FontWeight', 'bold');
ylim(ax5, [0 115]); grid(ax5, 'on');

% -- Panels 7-8 (bottom right, spanning 2 cols): TP/FP/FN grouped bar -----
ax6 = subplot(2, 4, [7 8]);
tp_fp_fn = zeros(3, 3);   % [method x TP/FP/FN]
for m = 1:3
    tp_fp_fn(m,1) = acc_data(m).tp;
    tp_fp_fn(m,2) = acc_data(m).fp;
    tp_fp_fn(m,3) = acc_data(m).fn;
end
b2 = bar(ax6, tp_fp_fn, 'grouped');
b2(1).FaceColor = [0.25 0.75 0.25];   % TP -- green
b2(2).FaceColor = [1.00 0.35 0.25];   % FP -- red
b2(3).FaceColor = [0.55 0.55 0.55];   % FN -- grey
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

% Save figure
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_path = fullfile(params.log_dir, ...
        sprintf('accuracy_histogram_%s.png', params.run_timestamp));
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Accuracy histogram figure saved -> accuracy_histogram_%s.png\n', ...
            params.run_timestamp);
end
end

%% ======================================================================
%  LOCAL HELPER FUNCTIONS -- Option 7 alignment pipeline (opt6 variants)
%  These are verbatim copies of the same-named functions in
%  run_segmentation_validation_10x.m, renamed with _opt6 suffix to avoid
%  MATLAB name collision when both files are on the path.
%% ======================================================================

function [mosaic, canvas_info, px_10x] = build_10x_mosaic_opt6(params)
data = bfopen(params.nd2_10x_path);
n_series = size(data, 1);
omeMeta  = data{1, 4};
try
    px_10x = double(omeMeta.getPixelsPhysicalSizeX(0).value());
catch
    px_10x = 0.655;
end
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
[H_tile, W_tile] = size(data{1,1}{params.channel_10x, 1});
stage_x_img = -stage_x;
stage_y_img = -stage_y;
stage_x_px = (stage_x_img - min(stage_x_img)) / px_10x;
stage_y_px = (stage_y_img - min(stage_y_img)) / px_10x;
canvas_W = round(max(stage_x_px) + W_tile);
canvas_H = round(max(stage_y_px) + H_tile);
mosaic   = zeros(canvas_H, canvas_W);
for s = 1:n_series
    raw = double(data{s,1}{params.channel_10x, 1});
    raw = raw ./ (max(raw(:)) + eps);
    tl_c = round(stage_x_px(s)) + 1;
    tl_r = round(stage_y_px(s)) + 1;
    r1 = max(1, tl_r);             r2 = min(canvas_H, tl_r + H_tile - 1);
    c1 = max(1, tl_c);             c2 = min(canvas_W, tl_c + W_tile - 1);
    mosaic(r1:r2, c1:c2) = max(mosaic(r1:r2, c1:c2), raw(1:r2-r1+1, 1:c2-c1+1));
end
canvas_info.stage_x_px  = stage_x_px;
canvas_info.stage_y_px  = stage_y_px;
canvas_info.stage_x_um  = stage_x;
canvas_info.stage_y_um  = stage_y;
canvas_info.tile_H      = H_tile;
canvas_info.tile_W      = W_tile;
canvas_info.canvas_H    = canvas_H;
canvas_info.canvas_W    = canvas_W;
canvas_info.x_origin_um = min(stage_x_img);
canvas_info.y_origin_um = min(stage_y_img);
fprintf('[10X]  %d FOVs placed on %dx%d canvas\n', n_series, canvas_W, canvas_H);
end


function [mips, coords, px_40x] = load_40x_mips_opt6(params)
data = bfopen(params.nd2_40x_path);
n_series_total = size(data, 1);
omeMeta = data{1, 4};
try
    px_40x = double(omeMeta.getPixelsPhysicalSizeX(0).value());
catch
    px_40x = 0.1625;
end
zpf   = params.z_planes_per_40x_fov;
n_fov = floor(n_series_total / zpf);
[H40, W40] = size(data{1,1}{params.channel_40x, 1});
mips   = cell(1, n_fov);
coords = struct('fov', num2cell(1:n_fov), ...
                'stage_x_um', num2cell(zeros(1,n_fov)), ...
                'stage_y_um', num2cell(zeros(1,n_fov)));
for f = 1:n_fov
    s_start = (f-1)*zpf + 1;
    s_end   = s_start + zpf - 1;
    mip = zeros(H40, W40);
    for s = s_start:s_end
        mip = max(mip, double(data{s,1}{params.channel_40x, 1}));
    end
    mips{f} = mip ./ (max(mip(:)) + eps);
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


function [gt_canvas, tl_c, tl_r] = place_gt_in_canvas_opt6(params, fov_40x_coords, canvas_info, px_10x, px_40x)
% Returns gt_canvas plus the top-left placement of the 40X mask inside it
% (tl_c, tl_r in canvas pixel coords).  Those are needed to map a 10X ROI
% back to the original 40X gt_mask coordinate frame.
gt_mask_raw  = uint16(imread(params.gt_mask_path, 1));
scale_factor = px_40x / px_10x;
gt_scaled    = imresize(gt_mask_raw, scale_factor, 'nearest');
[H_gs, W_gs] = size(gt_scaled);
fov1_stage_x_img = -fov_40x_coords(1).stage_x_um;
fov1_stage_y_img = -fov_40x_coords(1).stage_y_um;
fov1_cx_px = (fov1_stage_x_img - canvas_info.x_origin_um) / px_10x;
fov1_cy_px = (fov1_stage_y_img - canvas_info.y_origin_um) / px_10x;
tl_c = round(fov1_cx_px - W_gs/2 + canvas_info.tile_W/2) + 1;
tl_r = round(fov1_cy_px - H_gs/2 + canvas_info.tile_H/2) + 1;
gt_canvas = zeros(canvas_info.canvas_H, canvas_info.canvas_W, 'uint16');
r1 = max(1, tl_r); r2 = min(canvas_info.canvas_H, tl_r + H_gs - 1);
c1 = max(1, tl_c); c2 = min(canvas_info.canvas_W, tl_c + W_gs - 1);
gt_canvas(r1:r2, c1:c2) = gt_scaled((r1-tl_r+1):(r2-tl_r+1), (c1-tl_c+1):(c2-tl_c+1));
fprintf('[GT]   GT mask placed at canvas (%d,%d)\n', tl_c, tl_r);
end


function [crop_img, crop_mask, bbox] = crop_overlap_opt6(mosaic, gt_canvas)
[r_mask, c_mask] = find(gt_canvas > 0);
if isempty(r_mask)
    error('No GT mask overlap with 10X canvas. Check stage coordinates.');
end
r1 = min(r_mask); r2 = max(r_mask);
c1 = min(c_mask); c2 = max(c_mask);
[r_img, c_img] = find(mosaic > 0);
r1 = max(r1, min(r_img)); r2 = min(r2, max(r_img));
c1 = max(c1, min(c_img)); c2 = min(c2, max(c_img));
crop_img  = mosaic(r1:r2, c1:c2);
crop_mask = gt_canvas(r1:r2, c1:c2);
bbox      = [c1, r1, c2-c1+1, r2-r1+1];
end


function [shift_r, shift_c, peak_ncc] = refine_alignment_ncc_opt6(seg_bin, gt_bin, search_px)
[H, W] = size(seg_bin);
gt_pad = padarray(double(gt_bin), [search_px, search_px], 0, 'both');
C      = normxcorr2(double(seg_bin), gt_pad);
zr = search_px + H;
zc = search_px + W;
mask_ncc = false(size(C));
r_lo = max(1, zr - search_px); r_hi = min(size(C,1), zr + search_px);
c_lo = max(1, zc - search_px); c_hi = min(size(C,2), zc + search_px);
mask_ncc(r_lo:r_hi, c_lo:c_hi) = true;
Cm = C; Cm(~mask_ncc) = -inf;
[peak_ncc, idx] = max(Cm(:));
[pr, pc] = ind2sub(size(C), idx);
shift_r = zr - pr;
shift_c = zc - pc;
end


function plot_overlap_diagnostic_opt6(crop_10x, crop_gt, px_10x, params, tag)
if nargin < 5 || isempty(tag); tag = 'initial'; end
fig = figure('Name', ['Opt6 Overlap Diagnostic - ' tag], 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.75]);
sgtitle(sprintf('10X Crop vs 40X GT Alignment (%s)', tag), 'FontSize', 13, 'FontWeight', 'bold');
ax1 = subplot(1,3,1);
imagesc(ax1, crop_10x); colormap(ax1,'gray'); axis(ax1,'image'); axis(ax1,'off');
title(ax1, sprintf('10X crop (%dx%d px)', size(crop_10x,2), size(crop_10x,1)), 'FontSize', 10);
ax2 = subplot(1,3,2);
imagesc(ax2, crop_gt); colormap(ax2,'jet'); axis(ax2,'image'); axis(ax2,'off');
title(ax2, sprintf('40X GT mask in 10X scale (%d nuclei)', ...
      numel(unique(crop_gt(crop_gt>0)))), 'FontSize', 10);
ax3 = subplot(1,3,3);
rgb = repmat(mat2gray(crop_10x), [1 1 3]);
gt_edge = edge(crop_gt > 0, 'canny');
rgb(:,:,1) = min(1, rgb(:,:,1) + gt_edge);
imagesc(ax3, rgb); axis(ax3,'image'); axis(ax3,'off');
title(ax3, 'Overlay: 10X grayscale + GT boundaries (red)', 'FontSize', 10);
drawnow;
if isfield(params,'log_dir') && ~isempty(params.log_dir)
    exportgraphics(fig, fullfile(params.log_dir, ...
        sprintf('overlap_diagnostic_%s_%s.png', tag, params.run_timestamp)), 'Resolution', 150);
end
end
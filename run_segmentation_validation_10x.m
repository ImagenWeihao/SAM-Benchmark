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

if nargin < 1; error('params struct required'); end

fprintf('\n====================================================\n');
fprintf('  SEGMENTATION VALIDATION  10X input vs 40X GT mask\n');
fprintf('====================================================\n\n');

%% ── Step 1: Load 10X and build stitched mosaic ──────────────────────────
fprintf('[10X]  Loading overview and building mosaic...\n');
[mosaic_10x, canvas_info, px_10x] = build_10x_mosaic(params);
fprintf('[10X]  Mosaic: %d x %d px  (%.0f x %.0f um)\n', ...
        size(mosaic_10x,2), size(mosaic_10x,1), ...
        size(mosaic_10x,2)*px_10x, size(mosaic_10x,1)*px_10x);

%% ── Step 2: Load 40X FOVs, get stage coords, compute MIPs ───────────────
fprintf('\n[40X]  Loading FOVs + computing MIPs + reading stage coords...\n');
[fov_40x_mips, fov_40x_coords, px_40x] = load_40x_mips(params);

%% ── Step 3: Place 40X GT mask into 10X canvas ───────────────────────────
fprintf('\n[GT]   Placing 40X GT mask into 10X canvas...\n');
gt_mask_canvas = place_gt_in_canvas(params, fov_40x_coords, canvas_info, px_10x, px_40x);

%% ── Step 4: Crop overlap region ─────────────────────────────────────────
fprintf('\n[CROP] Cropping overlap region...\n');
[crop_10x, crop_gt, crop_bbox] = crop_overlap(mosaic_10x, gt_mask_canvas);
fprintf('[CROP] Overlap region: %d x %d px (%.0f x %.0f um)\n', ...
        size(crop_10x,2), size(crop_10x,1), ...
        size(crop_10x,2)*px_10x, size(crop_10x,1)*px_10x);
fprintf('[CROP] GT nuclei in overlap: %d\n', numel(unique(crop_gt(crop_gt>0))));

%% ── Step 5: Intermediate diagnostic (pre-refinement) ────────────────────
plot_overlap_diagnostic(crop_10x, crop_gt, px_10x, params);

%% ── Step 5.5: Large-template NCC refinement ─────────────────────────────
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

%% ── Step 6: Run segmentation methods on 10X overlap ─────────────────────
fprintf('\n[SEG]  Running segmentation on 10X overlap...\n');

% Method 1: MATLAB adaptive (reuse initial; identical to what we'd compute)
fprintf('[SEG]  Method 1: MATLAB adaptive segmentation (L1/L2)...\n');
seg_matlab = seg_initial;
stats_matlab = regionprops(seg_matlab, 'Area', 'Perimeter');
cells_matlab = struct('circularity', num2cell( ...
    arrayfun(@(s) 4*pi*s.Area / max(s.Perimeter^2, eps), stats_matlab)'));
fprintf('[SEG]    Detected %d nuclei\n', numel(cells_matlab));

% Method 2: Circularity-filtered (L2)
circ_thresh = 0.60;
fprintf('[SEG]  Method 2: Circularity-filtered (L2, thresh=%.2f)...\n', circ_thresh);
circ_vals    = [cells_matlab.circularity];
pass         = circ_vals >= circ_thresh;
seg_L2       = zeros(size(seg_matlab), 'uint16');
cids         = find(pass);
for k = 1:numel(cids); seg_L2(seg_matlab == cids(k)) = k; end
fprintf('[SEG]    %d/%d cells pass circularity\n', sum(pass), numel(cells_matlab));

% Method 3: StarDist
fprintf('[SEG]  Method 3: StarDist (L3)...\n');
seg_stardist = run_stardist_on_image(crop_10x);
fprintf('[SEG]    StarDist detected %d nuclei\n', max(seg_stardist(:)));

%% ── Step 7: Compute ROC vs 40X GT ───────────────────────────────────────
fprintf('\n[ROC]  Computing IoU-based ROC vs 40X GT mask...\n');
methods = {'MATLAB (L1/L2)', 'Circularity-filtered (L2)', 'StarDist (L3)'};
segs    = {seg_matlab, seg_L2, seg_stardist};
colors  = {[0.3 0.6 1.0], [1.0 0.4 0.3], [0.2 0.8 0.3]};
roc     = cell(1,3);

for m = 1:3
    roc{m} = compute_roc_objectlevel(crop_gt, segs{m});
    fprintf('[ROC]  %s: AUC=%.3f  F1@0.5=%.3f  detected=%d\n', ...
            methods{m}, roc{m}.auc, roc{m}.f1_at_half, roc{m}.n_pred);
end

%% ── Step 8: Display results ──────────────────────────────────────────────
plot_validation_results(crop_10x, crop_gt, segs, roc, methods, colors, px_10x, params);

%% ── Package output ───────────────────────────────────────────────────────
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
% Each FOV centred at (stage_x_s, stage_y_s) in µm
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
function gt_canvas = place_gt_in_canvas(params, fov_40x_coords, canvas_info, px_10x, px_40x)
%PLACE_GT_IN_CANVAS  Place 40X GT mask (FOV1) into the 10X stitched canvas.

% Load GT mask (Frame 0 = nuclei)
gt_mask = uint16(imread(params.gt_mask_path, 1));
[H_gt, W_gt] = size(gt_mask);

% Scale GT mask to 10X pixel size (nearest-neighbour to preserve labels)
scale_factor = px_40x / px_10x;
gt_scaled    = imresize(gt_mask, scale_factor, 'nearest');
[H_gs, W_gs] = size(gt_scaled);

% Place at FOV1 location in canvas (apply same transform as 10X)
fov1_stage_x_img = -fov_40x_coords(1).stage_x_um;
fov1_stage_y_img = -fov_40x_coords(1).stage_y_um;

% Convert to canvas pixel coords using canvas origin
fov1_cx_px = (fov1_stage_x_img - canvas_info.x_origin_um) / px_10x;
fov1_cy_px = (fov1_stage_y_img - canvas_info.y_origin_um) / px_10x;

% Top-left corner for placement (centred)
tl_c = round(fov1_cx_px - W_gs/2 + canvas_info.tile_W/2) + 1;
tl_r = round(fov1_cy_px - H_gs/2 + canvas_info.tile_H/2) + 1;

gt_canvas = zeros(canvas_info.canvas_H, canvas_info.canvas_W, 'uint16');
r1 = max(1, tl_r); r2 = min(canvas_info.canvas_H, tl_r + H_gs - 1);
c1 = max(1, tl_c); c2 = min(canvas_info.canvas_W, tl_c + W_gs - 1);
src_r = (r1-tl_r+1):(r2-tl_r+1);
src_c = (c1-tl_c+1):(c2-tl_c+1);

gt_canvas(r1:r2, c1:c2) = gt_scaled(src_r, src_c);
fprintf('[GT]   GT mask (%dx%d) placed at canvas (%d,%d)\n', ...
        W_gs, H_gs, tl_c, tl_r);
end


%% ======================================================================
function [crop_img, crop_mask, bbox] = crop_overlap(mosaic, gt_canvas)
%CROP_OVERLAP  Extract region where both 10X image and GT mask exist.
[r_mask, c_mask] = find(gt_canvas > 0);
if isempty(r_mask)
    error('No GT mask overlap with 10X canvas. Check stage coordinates.');
end
r1 = min(r_mask); r2 = max(r_mask);
c1 = min(c_mask); c2 = max(c_mask);

% Also clip to mosaic non-zero region
[r_img, c_img] = find(mosaic > 0);
r1 = max(r1, min(r_img)); r2 = min(r2, max(r_img));
c1 = max(c1, min(c_img)); c2 = min(c2, max(c_img));

crop_img  = mosaic(r1:r2, c1:c2);
crop_mask = gt_canvas(r1:r2, c1:c2);
bbox      = [c1, r1, c2-c1+1, r2-r1+1];
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
function plot_validation_results(crop_10x, crop_gt, segs, roc, methods, colors, px_10x, params)
%PLOT_VALIDATION_RESULTS  Segmentation overlays + ROC + morphology bars.
fig = figure('Name','10X vs 40X GT Validation','NumberTitle','off', ...
             'Units','normalized','Position',[0.02 0.02 0.96 0.92]);
sgtitle('Segmentation Validation  10X input vs 40X GT mask', ...
        'FontSize',13,'FontWeight','bold');

method_short = {'MATLAB (L1/L2)', 'Circ.-filter (L2)', 'StarDist (L3)'};

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
    n_missed = roc{m}.n_gt - round(roc{m}.tpr(round(end/2)) * roc{m}.n_gt);
    title(ax, sprintf('%s\nAUC=%.3f  F1=%.3f  detected=%d', ...
          method_short{m}, roc{m}.auc, roc{m}.f1_at_half, roc{m}.n_pred), ...
          'FontSize', 9);
end

% Panel 4: ROC curves
ax4 = subplot(2,3,4);
hold(ax4,'on');
plot(ax4, [0 1], [0 1], 'k--', 'DisplayName','Random');
for m = 1:3
    [fs, idx] = sort(roc{m}.fpr);
    ts = roc{m}.tpr(idx);
    plot(ax4, fs, ts, '-', 'Color', colors{m}, 'LineWidth', 2, ...
         'DisplayName', sprintf('%s (AUC=%.3f)', method_short{m}, roc{m}.auc));
end
xlabel(ax4, 'FPR (object-level)'); ylabel(ax4, 'TPR (Recall)');
title(ax4, 'ROC Curves', 'FontSize', 10);
legend(ax4, 'Location', 'southeast','FontSize',8);
xlim(ax4,[0 1]); ylim(ax4,[0 1]); grid(ax4,'on'); axis(ax4,'square');

% Panel 5: Detection summary
ax5 = subplot(2,3,5);
bar_data = zeros(3, 3);
for m = 1:3
    bar_data(m,1) = roc{m}.n_gt;
    bar_data(m,2) = roc{m}.n_pred;
    bar_data(m,3) = round(roc{m}.f1_at_half * roc{m}.n_gt);  % approx TP
end
b = bar(ax5, bar_data');
for m = 1:3; b(m).FaceColor = colors{m}; end
set(ax5, 'XTickLabel', {'GT nuclei','Predicted','TP at IoU=0.5'}, 'XTick', 1:3);
legend(ax5, method_short, 'Location', 'best', 'FontSize', 8);
title(ax5, 'Detection Counts', 'FontSize', 10);
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
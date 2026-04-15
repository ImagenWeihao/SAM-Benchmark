function results = run_nucleus_normality_test(img_data, params)
%RUN_NUCLEUS_NORMALITY_TEST  Two-stage nucleus normality classification pipeline.
%
%  Stage 1: Detect and segment all nuclei from DAPI (405nm) channel
%           Uses segment_nuclei.m (existing shared utility)
%  Stage 2: Classify each nucleus crop using lightweight UNet:
%           - normal          : round, compact, single nucleus
%           - abnormal_shape  : multi-lobular, blebbing (classes B, F)
%           - abnormal_count  : binucleated, polyploid, micronuclei (C, D, E)
%
%  Can run in two modes set by params.normality_mode:
%    'mock'   — rule-based proxy using circularity + area (no model needed)
%    'unet'   — uses trained/untrained UNet via Deep Learning Toolbox
%
%  Usage:
%    results = run_nucleus_normality_test(img_data, params)
%
%  Reference morphologies (from attached diagram):
%    A - Normal cell:          single round nucleus, 2N
%    B - Multi-lobular:        irregular/lobular nucleus boundary
%    C - Binucleated:          two separate nuclei of similar size
%    D - Polyploid:            single oversized nucleus (4N+)
%    E - Micronuclei:          main nucleus + small satellite nucleus
%    F - Blebbing nucleus:     protrusion/bleb on nucleus boundary
%    G - (advanced) nuclear substructure — requires higher resolution

fprintf('\n==============================================\n');
fprintf('  NUCLEUS NORMALITY CLASSIFICATION TEST\n');
fprintf('  Mode:    %s\n', params.normality_mode);
fprintf('  Channel: 405nm (DAPI)\n');
fprintf('==============================================\n\n');
t_start = tic;

%% ── Class definitions ────────────────────────────────────────────────────
CLASS_NORMAL   = 1;   % A
CLASS_ABN_SHAPE= 2;   % B, F
CLASS_ABN_COUNT= 3;   % C, D, E

class_names = {'normal', 'abnormal_shape', 'abnormal_count'};
class_desc  = { ...
    'Round, compact, single nucleus (A)', ...
    'Multi-lobular or blebbing nucleus (B, F)', ...
    'Binucleated, polyploid, or micronuclei (C, D, E)'};

%% ── Step 1: Segment all nuclei from DAPI channel ─────────────────────────
fprintf('[NNT] Step 1/4 — Segmenting nuclei from 405nm channel...\n');
nucleus_img = img_data.nucleus;   % already extracted as channel 1

% Use tighter params for normality testing (detect micronuclei too)
seg_params                      = params;
seg_params.min_nucleus_area     = 50;    % smaller minimum to catch micronuclei
[all_cells, seg_mask, bw_nuclei]= segment_nuclei_normality(nucleus_img, seg_params);

n_detected = numel(all_cells);
fprintf('[NNT] Detected %d nucleus objects\n', n_detected);
if n_detected == 0
    error('No nuclei detected. Check image and segmentation parameters.');
end

%% ── Step 2: Extract per-nucleus feature crops ────────────────────────────
fprintf('[NNT] Step 2/4 — Extracting nucleus crops (64x64 px)...\n');
crop_size = 64;
[H, W]    = size(nucleus_img);
half      = crop_size / 2;

crops = zeros(crop_size, crop_size, 1, n_detected, 'single');
for k = 1:n_detected
    cx = round(all_cells(k).centroid(1));
    cy = round(all_cells(k).centroid(2));
    x1 = max(1, cx-half+1); x2 = min(W, cx+half);
    y1 = max(1, cy-half+1); y2 = min(H, cy+half);
    crop = nucleus_img(y1:y2, x1:x2);

    % Pad to exact crop_size if at image edge
    padded = zeros(crop_size, crop_size, 'single');
    padded(1:size(crop,1), 1:size(crop,2)) = single(crop);
    crops(:,:,1,k) = padded;
end
fprintf('[NNT] Crops extracted: %dx%dx1x%d\n', crop_size, crop_size, n_detected);

%% ── Step 3: Classify each nucleus ───────────────────────────────────────
fprintf('[NNT] Step 3/4 — Classifying nuclei (mode: %s)...\n', params.normality_mode);

switch lower(params.normality_mode)
    case {'stardist', 'cellpose', 'nuclear_seg'}
        unet_params = params;
        unet_params.unet_model  = params.normality_mode;
        unet_params.unet_target = 'normal_nucleus';
        [classified, ~] = classify_cells_unet(all_cells, img_data, unet_params);
        class_ids     = ones(n_detected, 1, 'int32');
        confidences   = zeros(n_detected, 1);
        class_reasons = cell(n_detected, 1);
        for k = 1:n_detected
            switch classified(k).unet_class
                case 'abnormal_shape'; class_ids(k) = int32(2);
                case 'abnormal_count'; class_ids(k) = int32(3);
                otherwise;             class_ids(k) = int32(1);
            end
            confidences(k)   = classified(k).unet_score;
            class_reasons{k} = sprintf('%s: score=%.3f', params.normality_mode, classified(k).unet_score);
        end

    case 'mock'
        [class_ids, confidences, class_reasons] = ...
            classify_mock(all_cells, crops, params);

    case 'unet'
        [class_ids, confidences, class_reasons] = ...
            classify_unet(all_cells, crops, params);

    otherwise
        error('Unknown normality_mode: "%s". Use ''mock'', ''unet'', ''stardist'', ''cellpose'', or ''nuclear_seg''.', ...
              params.normality_mode);
end

%% ── Step 4: Package and report results ──────────────────────────────────
fprintf('\n[NNT] Step 4/4 — Compiling results...\n\n');

% Attach classification to cell structs
for k = 1:n_detected
    all_cells(k).norm_class      = class_ids(k);
    all_cells(k).norm_class_name = class_names{class_ids(k)};
    all_cells(k).norm_confidence = confidences(k);
    all_cells(k).norm_reason     = class_reasons{k};
end

% Count per class
n_normal     = sum(class_ids == CLASS_NORMAL);
n_abn_shape  = sum(class_ids == CLASS_ABN_SHAPE);
n_abn_count  = sum(class_ids == CLASS_ABN_COUNT);
n_abnormal   = n_abn_shape + n_abn_count;
pct_normal   = 100 * n_normal / n_detected;

% Print classification table
fprintf('%-6s  %-16s  %-8s  %-8s  %s\n', ...
        'CellID', 'Class', 'Confid.', 'Circ.', 'Reason');
fprintf('%s\n', repmat('-', 1, 72));
for k = 1:n_detected
    c = all_cells(k);
    marker = '  ';
    if c.norm_class ~= CLASS_NORMAL; marker = '**'; end
    fprintf('%s%-4d  %-16s  %-8.3f  %-8.3f  %s\n', ...
            marker, c.id, c.norm_class_name, c.norm_confidence, ...
            c.circularity, c.norm_reason);
end

fprintf('\n[NNT] -- Classification Summary ------------------\n');
fprintf('[NNT]   Total nuclei:     %d\n',         n_detected);
fprintf('[NNT]   Normal:           %d (%.1f%%)\n', n_normal,    pct_normal);
fprintf('[NNT]   Abnormal shape:   %d (%.1f%%)\n', n_abn_shape, 100*n_abn_shape/n_detected);
fprintf('[NNT]   Abnormal count:   %d (%.1f%%)\n', n_abn_count, 100*n_abn_count/n_detected);
fprintf('[NNT]   Total abnormal:   %d (%.1f%%)\n', n_abnormal,  100*n_abnormal/n_detected);
fprintf('[NNT] -------------------------------------------------\n\n');

elapsed = toc(t_start);

%% ── Package output struct ────────────────────────────────────────────────
results.cells         = all_cells;
results.seg_mask      = seg_mask;
results.bw_nuclei     = bw_nuclei;
results.crops         = crops;
results.class_ids     = class_ids;
results.confidences   = confidences;
results.class_names   = class_names;
results.class_desc    = class_desc;
results.n_detected    = n_detected;
results.n_normal      = n_normal;
results.n_abn_shape   = n_abn_shape;
results.n_abn_count   = n_abn_count;
results.pct_normal    = pct_normal;
results.elapsed_sec   = elapsed;
results.mode          = params.normality_mode;

%% ── Visualise ────────────────────────────────────────────────────────────
plot_normality_results(results, nucleus_img, params);
end


%% ======================================================================
%  CLASSIFIER: MOCK (rule-based, no model needed)
%  Uses circularity, area, and solidity to mimic the 3 classes
%% ======================================================================
function [class_ids, confidences, reasons] = classify_mock(cells, crops, params)
%CLASSIFY_MOCK  Rule-based proxy classifier using DAPI morphology features.
%
%  Rules derived from the reference morphology diagram:
%
%  abnormal_count (C/D/E):
%    - Very large nucleus (area > 2.5x mean): likely polyploid (D)
%    - Very small nucleus (area < 0.3x mean): likely micronucleus (E)
%    - Low solidity with medium area: may indicate binucleation (C)
%      [True binucleation detection requires proximity analysis of
%       nearby same-size nuclei — flagged here via area clustering]
%
%  abnormal_shape (B/F):
%    - Low circularity (< threshold): multi-lobular or blebbing
%    - Low solidity with near-average area: blebbing / protrusion
%
%  normal (A):
%    - High circularity AND solidity AND near-average area

n       = numel(cells);
areas   = [cells.area];
mean_a  = mean(areas);

class_ids   = ones(n, 1, 'int32');   % default: normal
confidences = zeros(n, 1);
reasons     = cell(n, 1);

circ_thresh  = params.normality_circ_threshold;    % default 0.80
solid_thresh = params.normality_solid_threshold;   % default 0.85
area_hi      = params.normality_area_ratio_hi;     % default 2.5
area_lo      = params.normality_area_ratio_lo;     % default 0.30

for k = 1:n
    c     = cells(k);
    circ  = c.circularity;
    solid = c.solidity;
    area  = c.area;
    a_rat = area / mean_a;

    if a_rat > area_hi
        % Polyploid (D): nucleus much larger than average
        class_ids(k)   = int32(3);   % abnormal_count
        confidences(k) = min(1.0, 0.6 + (a_rat - area_hi) * 0.1);
        reasons{k}     = sprintf('Oversized nucleus (%.1fx mean area) — possible polyploid (D)', a_rat);

    elseif a_rat < area_lo
        % Micronucleus (E): very small object
        class_ids(k)   = int32(3);   % abnormal_count
        confidences(k) = min(1.0, 0.7 + (area_lo - a_rat) * 0.5);
        reasons{k}     = sprintf('Undersized nucleus (%.1fx mean area) — possible micronucleus (E)', a_rat);

    elseif solid < solid_thresh && a_rat > 0.6 && a_rat < 1.8 && circ > circ_thresh
        % Binucleated (C): moderate area, reasonable circularity but low solidity
        % (two nuclei touching → concave hull)
        class_ids(k)   = int32(3);   % abnormal_count
        confidences(k) = 0.5 + (solid_thresh - solid) * 1.5;
        confidences(k) = min(1.0, confidences(k));
        reasons{k}     = sprintf('Low solidity (%.2f) with normal area — possible binucleation (C)', solid);

    elseif circ < circ_thresh
        % Multi-lobular / blebbing (B/F)
        class_ids(k)   = int32(2);   % abnormal_shape
        confidences(k) = min(1.0, 0.55 + (circ_thresh - circ) * 2.0);
        if circ < circ_thresh - 0.15
            reasons{k} = sprintf('Very low circularity (%.2f) — likely multi-lobular nucleus (B)', circ);
        else
            reasons{k} = sprintf('Low circularity (%.2f) — possible blebbing nucleus (F)', circ);
        end

    elseif solid < solid_thresh && circ < circ_thresh + 0.05
        % Blebbing (F): near-normal circularity but low solidity
        class_ids(k)   = int32(2);   % abnormal_shape
        confidences(k) = min(1.0, 0.5 + (solid_thresh - solid) * 2.0);
        reasons{k}     = sprintf('Low solidity (%.2f) with marginal circularity — possible blebbing (F)', solid);

    else
        % Normal (A)
        class_ids(k)   = int32(1);
        confidences(k) = (circ + solid) / 2 * min(1.0, 1 - abs(a_rat - 1.0) * 0.3);
        confidences(k) = max(0.5, min(1.0, confidences(k)));
        reasons{k}     = sprintf('Normal: circ=%.2f, solid=%.2f, area_ratio=%.2f', circ, solid, a_rat);
    end
end
end


%% ======================================================================
%  CLASSIFIER: UNET (Deep Learning Toolbox)
%% ======================================================================
function [class_ids, confidences, reasons] = classify_unet(cells, crops, params)
%CLASSIFY_UNET  Run lightweight UNet on nucleus crops.
%
%  If a trained model exists at params.normality_model_path, it is loaded.
%  Otherwise builds the untrained skeleton and falls back to mock scoring
%  with a warning — this ensures the pipeline always produces output.

n         = numel(cells);
class_ids = ones(n, 1, 'int32');
confidences = zeros(n, 1);
reasons   = cell(n, 1);

% Try to load trained model
model_loaded = false;
if isfield(params, 'normality_model_path') && ~isempty(params.normality_model_path) ...
        && isfile(params.normality_model_path)
    fprintf('[NNT]  Loading trained UNet from: %s\n', params.normality_model_path);
    loaded = load(params.normality_model_path, 'net');
    net    = loaded.net;
    model_loaded = true;
else
    fprintf('[NNT]  No trained model found. Building architecture skeleton...\n');
    net = build_nucleus_classifier_unet([64 64 1], 3);
    fprintf('[NNT]  WARNING: Model is untrained — falling back to mock scoring.\n');
    fprintf('[NNT]  To train: call train_nucleus_classifier() and save the net.\n');
end

if model_loaded
    % Run semantic segmentation on all crops
    for k = 1:n
        crop = crops(:,:,:,k);
        try
            seg_map     = semanticseg(crop, net);
            % Dominant class in centre region (nucleus body)
            ctr         = round(size(seg_map,1)/2);
            region      = seg_map(ctr-8:ctr+8, ctr-8:ctr+8);
            counts      = histcounts(uint8(region), 1:4);
            [conf, cls] = max(counts);
            class_ids(k)   = int32(cls);
            confidences(k) = conf / sum(counts);
            reasons{k}     = sprintf('UNet segmentation: dominant class %d in centre region', cls);
        catch
            class_ids(k)   = int32(1);
            confidences(k) = 0.5;
            reasons{k}     = 'UNet inference failed — defaulting to normal';
        end
    end
else
    % Fall back to mock scoring
    [class_ids, confidences, reasons] = classify_mock(cells, crops, params);
    reasons = cellfun(@(r) ['[skeleton fallback] ' r], reasons, 'UniformOutput', false);
end
end


%% ======================================================================
%  SPECIALISED SEGMENTATION for normality testing
%  (lower min area than standard to catch micronuclei)
%% ======================================================================
function [cells, seg_mask, bw_nuclei] = segment_nuclei_normality(nucleus_img, params)
%SEGMENT_NUCLEI_NORMALITY  Modified segmentation tuned for normality testing.
%  Extends segment_nuclei.m with:
%    - Smaller minimum area (catches micronuclei, class E)
%    - Two-pass detection: first normal nuclei, then small satellites

% Primary segmentation (standard)
[cells, seg_mask, bw_nuclei] = segment_nuclei(nucleus_img, params);

% Second pass: detect micronuclei (objects too small for primary pass)
% Threshold: objects between 50-200 px² that weren't already labelled
img_smooth = imgaussfilt(nucleus_img, 1.0);
bw_all     = imbinarize(img_smooth, 'adaptive', 'Sensitivity', 0.45, ...
                        'ForegroundPolarity', 'bright');
bw_all     = imfill(bw_all, 'holes');

% Remove already-detected regions
bw_new = bw_all & ~(seg_mask > 0);
bw_new = bwareaopen(bw_new, 50);    % min 50 px²

% Remove anything larger than a micronucleus threshold
micro_max = 200;
bw_micro  = bw_new & ~bwareaopen(bw_new, micro_max);

if any(bw_micro(:))
    [micro_label, n_micro] = bwlabel(bw_micro);
    micro_props = regionprops(micro_label, nucleus_img, ...
        'Area', 'Centroid', 'Perimeter', 'BoundingBox', 'MeanIntensity');

    n_existing = numel(cells);
    for m = 1:n_micro
        % Build new_cell using existing cell as template to match field structure
        new_cell             = cells(1);   % copy all fields from first cell
        new_cell.id          = n_existing + m;
        new_cell.centroid    = micro_props(m).Centroid;
        new_cell.area        = micro_props(m).Area;
        new_cell.perimeter   = micro_props(m).Perimeter;
        new_cell.bbox        = micro_props(m).BoundingBox;
        new_cell.mean_int    = micro_props(m).MeanIntensity;
        p = new_cell.perimeter;
        a = new_cell.area;
        new_cell.circularity = (p > 0) * min(1, (4*pi*a) / (p^2 + eps));
        cells(end+1) = new_cell; %#ok<AGROW>
    end

    % Update seg_mask
    offset            = max(seg_mask(:));
    micro_label(micro_label > 0) = micro_label(micro_label > 0) + double(offset);
    seg_mask          = seg_mask + uint16(micro_label);

    fprintf('[SEG]  Micronuclei pass: found %d additional small objects\n', n_micro);
end
end


%% ======================================================================
%  VISUALISATION
%% ======================================================================
function plot_normality_results(R, nucleus_img, params)

class_colors = [0.2 0.8 0.2;   % normal       — green
                1.0 0.5 0.0;   % abn_shape    — orange
                0.9 0.1 0.1];  % abn_count    — red

fig = figure('Name', 'Nucleus Normality Classification', 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.03 0.05 0.94 0.88]);
sgtitle(sprintf('Nucleus Normality Test — 405nm DAPI Channel  |  Mode: %s', R.mode), ...
        'FontSize', 15, 'FontWeight', 'bold');

%% Panel 1 -- Raw DAPI image with all nuclei overlaid
ax1 = subplot(2, 4, 1);
imagesc(ax1, nucleus_img); colormap(ax1, 'gray'); axis(ax1,'image'); axis(ax1,'off');
hold(ax1, 'on');
for k = 1:R.n_detected
    c     = R.cells(k).centroid;
    cls   = R.class_ids(k);
    col   = class_colors(cls,:);
    plot(ax1, c(1), c(2), 'o', 'Color', col, 'MarkerFaceColor', col, ...
         'MarkerSize', 10, 'LineWidth', 1.5);
    text(ax1, c(1)+6, c(2), num2str(R.cells(k).id), 'Color', 'white', 'FontSize', 7);
end
% Fix 1: full legend labels
title(ax1, {sprintf('DAPI — All nuclei (n=%d)', R.n_detected), ...
            'green=normal  orange=abnormal shape  red=abnormal count'}, 'FontSize', 8);

%% Panel 2 -- Segmentation mask coloured by class
ax2 = subplot(2, 4, 2);
overlay = zeros(size(nucleus_img,1), size(nucleus_img,2), 3);
for k = 1:R.n_detected
    % Fix 2: match seg_mask label to cell index (1-based ID)
    mask_k = (R.seg_mask == k);
    col    = class_colors(R.class_ids(k),:);
    for ch = 1:3
        overlay(:,:,ch) = overlay(:,:,ch) + mask_k * col(ch);
    end
end
% Dim background + coloured overlay
bg = repmat(nucleus_img, [1 1 3]) * 0.25;
imagesc(ax2, min(1, bg + overlay)); axis(ax2,'image'); axis(ax2,'off');
% Fix 2: add clear legend
legend_entries = {'Normal', 'Abn. Shape', 'Abn. Count'};
for cls = 1:3
    patch(ax2, NaN, NaN, class_colors(cls,:), 'DisplayName', legend_entries{cls});
end
legend(ax2, 'show', 'Location', 'southeast', 'FontSize', 7);
title(ax2, {'Segmentation', 'coloured by class'}, 'FontSize', 10);

%% Panel 3 -- Pie chart of classification breakdown
ax3 = subplot(2, 4, 3);
counts   = [R.n_normal, R.n_abn_shape, R.n_abn_count];
lbls     = {sprintf('Normal (n=%d)',       R.n_normal), ...
            sprintf('Abn. Shape (n=%d)',   R.n_abn_shape), ...
            sprintf('Abn. Count (n=%d)',   R.n_abn_count)};
non_zero = counts > 0;
if sum(non_zero) > 0
    pie(ax3, counts(non_zero), lbls(non_zero));
    colormap(ax3, class_colors(non_zero,:));
end
% Fix 3: clean two-line title, no overlap
title(ax3, {'Classification breakdown', sprintf('%.0f%% normal', R.pct_normal)}, 'FontSize', 10);

%% Panel 4 -- Circularity vs Solidity scatter, coloured by class
ax4 = subplot(2, 4, 4);
hold(ax4, 'on');
class_labels_leg = {'Normal (A)', 'Abn. Shape (B/F)', 'Abn. Count (C/D/E)'};
leg_h    = gobjects(3,1);
all_circ = [R.cells.circularity];
all_solid= ones(1, R.n_detected);   % default if solidity missing
if isfield(R.cells, 'solidity')
    all_solid = [R.cells.solidity];
end

for cls = 1:3
    idx = find(R.class_ids == cls);
    if ~isempty(idx)
        circ_v  = all_circ(idx);
        solid_v = all_solid(idx);
        leg_h(cls) = scatter(ax4, circ_v, solid_v, 60, class_colors(cls,:), ...
                             'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        for i = 1:numel(idx)
            text(ax4, circ_v(i)+0.003, solid_v(i), ...
                 num2str(R.cells(idx(i)).id), 'FontSize', 7, 'Color', [0.2 0.2 0.2]);
        end
    else
        leg_h(cls) = scatter(ax4, NaN, NaN, 1, class_colors(cls,:), 'filled');
    end
end

xline(ax4, params.normality_circ_threshold, '--k', 'LineWidth', 1, ...
      'Label', sprintf('circ=%.2f', params.normality_circ_threshold));
yline(ax4, params.normality_solid_threshold, ':k', 'LineWidth', 1, ...
      'Label', sprintf('solid=%.2f', params.normality_solid_threshold));
xlabel(ax4, 'Circularity'); ylabel(ax4, 'Solidity');
legend(ax4, leg_h, class_labels_leg, 'Location', 'best', 'FontSize', 8);
title(ax4, {'Feature Space', '(circularity vs solidity)'}, 'FontSize', 10);

% Fix 4: auto-zoom to data range with small padding
circ_pad  = 0.05;
solid_pad = 0.05;
xlim(ax4, [max(0.5, min(all_circ)-circ_pad),  min(1.02, max(all_circ)+circ_pad)]);
ylim(ax4, [max(0.3, min(all_solid)-solid_pad), min(1.05, max(all_solid)+solid_pad)]);
grid(ax4, 'on');

%% Panels 5-8 -- Crop gallery: show up to 8 cells, coloured border by class
n_show = min(R.n_detected, 4);
% Sort: show abnormal first
[~, sort_ord] = sort(R.class_ids(1:R.n_detected), 'descend');
show_idx = sort_ord(1:n_show);

for g = 1:n_show
    ax = subplot(2, 4, 4+g);
    k   = show_idx(g);
    cls = R.class_ids(k);
    col = class_colors(cls,:);

    crop_k = R.crops(:,:,1,k);
    imagesc(ax, crop_k); colormap(ax, 'gray'); axis(ax,'image'); axis(ax,'off');

    % Coloured border to indicate class
    set(ax, 'XColor', col, 'YColor', col, 'LineWidth', 3);
    axis(ax, 'on'); ax.XTick = []; ax.YTick = [];

    conf = R.confidences(k);
    % Use short class labels to avoid truncation
    switch R.cells(k).norm_class_name
        case 'normal';          cls_label = 'Normal';
        case 'abnormal_shape';  cls_label = 'Abn. Shape';
        case 'abnormal_count';  cls_label = 'Abn. Count';
        otherwise;              cls_label = R.cells(k).norm_class_name;
    end
    title(ax, {sprintf('Cell %d: %s', R.cells(k).id, cls_label), ...
               sprintf('conf=%.2f  circ=%.2f', conf, R.cells(k).circularity)}, ...
          'FontSize', 8, 'Color', col*0.7);
end

drawnow;

%% Save figure
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_path = fullfile(params.log_dir, ...
        sprintf('normality_test_%s.png', params.run_timestamp));
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Normality figure saved.\n');
end
end
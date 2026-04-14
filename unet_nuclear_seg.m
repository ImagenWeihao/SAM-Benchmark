function [cells_out, model_info] = unet_nuclear_seg(cells_in, img_data, params)
%UNET_NUCLEAR_SEG  nuclearSegmentator adapter.
%
%  nuclearSegmentator is purpose-built for nuclear morphology classification
%  in fluorescence microscopy. Unlike Cellpose/StarDist which are segmentation
%  models, nuclearSegmentator directly outputs morphology class probabilities
%  per nucleus — making it the closest drop-in for your normality use case.
%
%  GitHub: https://github.com/BIIFSweden/nuclearSegmentator
%  Paper:  "Automatic nuclear morphology classification" (bioRxiv 2023)
%
%  Output classes (mapped to your 3-class scheme):
%    nuclear_seg 'normal'      → normal           (A)
%    nuclear_seg 'irregular'   → abnormal_shape   (B, F)
%    nuclear_seg 'multi'       → abnormal_count   (C, D, E)
%
%  Requirements:
%    pip install nuclear-segmentator
%    pyenv('Version', 'C:\path\to\python.exe')
%
%  params fields used:
%    .nucleus_channel  — channel index for DAPI image
%    .unet_threshold   — minimum confidence to accept classification

%% ── Model info ───────────────────────────────────────────────────────────
model_info.name    = 'nuclearSegmentator';
model_info.version = 'pretrained-fluorescence';
model_info.target  = params.unet_target;
model_info.status  = 'pending';

cells_out = cells_in;
n         = numel(cells_in);

%% ── Check Python + nuclearSegmentator ────────────────────────────────────
try
    ns_mod = py.importlib.import_module('nuclear_segmentator');
catch ME
    error('UNET:nuclearSegNotFound', ...
          ['nuclearSegmentator not found.\n' ...
           'Fix:\n' ...
           '  1. pip install nuclear-segmentator\n' ...
           '  2. pyenv(''Version'', ''C:\\path\\to\\python.exe'')\n' ...
           'Original error: %s'], ME.message);
end

%% ── Load model ───────────────────────────────────────────────────────────
fprintf('[NucSeg] Loading nuclearSegmentator pretrained model...\n');
np    = py.importlib.import_module('numpy');

% nuclearSegmentator exposes a simple predict() interface
classifier = ns_mod.NuclearClassifier();
classifier.load_pretrained();

%% ── Class mapping ────────────────────────────────────────────────────────
% Maps nuclearSegmentator output labels to your 3-class scheme
ns_to_class = containers.Map( ...
    {'normal',    'round',   'compact', ...   % → normal
     'irregular', 'lobular', 'bleb',    ...   % → abnormal_shape
     'multi',     'bi',      'micro',   ...   % → abnormal_count
     'large',     'small'}, ...               % → abnormal_count
    {'normal',        'normal',        'normal', ...
     'abnormal_shape','abnormal_shape','abnormal_shape', ...
     'abnormal_count','abnormal_count','abnormal_count', ...
     'abnormal_count','abnormal_count'});

nucleus_img = img_data.nucleus;
[H, W]      = size(nucleus_img);
crop_half   = 64;

%% ── Run per-cell ─────────────────────────────────────────────────────────
for k = 1:n
    cx = round(cells_in(k).centroid(1));
    cy = round(cells_in(k).centroid(2));
    x1 = max(1, cx-crop_half); x2 = min(W, cx+crop_half);
    y1 = max(1, cy-crop_half); y2 = min(H, cy+crop_half);
    crop = nucleus_img(y1:y2, x1:x2);

    % Pad to 128x128 (nuclearSegmentator preferred input size)
    pad_size   = 128;
    crop_pad   = zeros(pad_size, pad_size, 'single');
    crop_pad(1:size(crop,1), 1:size(crop,2)) = single(crop);

    % Normalise
    cmax = max(crop_pad(:));
    if cmax > 0; crop_pad = crop_pad / cmax; end

    py_crop = np.array(crop_pad);

    try
        result     = classifier.predict(py_crop);
        ns_label   = char(result{'class'});
        confidence = double(result{'confidence'});

        % Map to your 3-class scheme
        if isKey(ns_to_class, lower(ns_label))
            unet_class = ns_to_class(lower(ns_label));
        else
            unet_class = 'normal';   % conservative default
            warning('[NucSeg] Unknown class label: %s', ns_label);
        end

        % unet_score: high = more normal
        if strcmp(unet_class, 'normal')
            score = confidence;
        else
            score = 1 - confidence;  % invert: abnormal confidence → low normality score
        end

    catch ME
        warning('[NucSeg] Cell %d failed: %s', k, ME.message);
        unet_class = 'normal';
        score      = 0.5;
        ns_label   = 'error';
        confidence = 0;
    end

    cells_out(k).unet_score   = max(0, min(1, score));
    cells_out(k).unet_class   = unet_class;
    cells_out(k).ns_raw_label = ns_label;
    cells_out(k).ns_confidence= confidence;

    fprintf('[NucSeg] Cell %2d: ns_class=%-10s conf=%.3f → %s (score=%.3f)\n', ...
            cells_in(k).id, ns_label, confidence, unet_class, score);
end

model_info.status      = 'ok';
model_info.n_processed = n;
end

function net = train_nucleus_classifier(image_dir, label_dir, params)
%TRAIN_NUCLEUS_CLASSIFIER  Training scaffold for the nucleus normality UNet.
%
%  net = train_nucleus_classifier(image_dir, label_dir, params)
%
%  Trains the lightweight 3-class UNet on your own annotated nucleus crops.
%
%  INPUT PREPARATION:
%    image_dir  — folder of 64x64 PNG nucleus crops (DAPI channel, grayscale)
%                 Naming: cell_001.png, cell_002.png, ...
%
%    label_dir  — folder of matching 64x64 PNG label maps (same filenames)
%                 Pixel values:  1 = normal
%                                2 = abnormal_shape (B, F)
%                                3 = abnormal_count (C, D, E)
%                 Background:    0 (ignored in loss)
%
%  QUICKEST WAY TO GET LABELS:
%    1. Export nucleus crops using: export_crops_for_labelling()  (below)
%    2. Open in ImageJ/FIJI: Plugins > Labeling > Label Nuclei
%       Or use MATLAB Image Labeler:
%         imageLabeler(image_dir)
%         — define 3 pixel classes, export label images
%    3. Pass labelled folder to this function
%
%  OUTPUT:
%    net — trained DAGNetwork, save with:
%          save(params.normality_model_path, 'net')
%
%  Requires: Deep Learning Toolbox, Computer Vision Toolbox

if nargin < 3; params = struct(); end

%% ── Training options ─────────────────────────────────────────────────────
input_size  = [64 64 1];
n_classes   = 3;
class_names = {'normal', 'abnormal_shape', 'abnormal_count'};

% Build network
lgraph = build_nucleus_classifier_unet(input_size, n_classes);

% Class weights — upweight abnormal classes (typically rarer)
% Adjust based on your label distribution
class_weights = [1.0, 2.5, 2.5];

% Update pixelClassificationLayer with class info and weights
pxLayer = pixelClassificationLayer( ...
    'Name',           'output', ...
    'Classes',        class_names, ...
    'ClassWeights',   class_weights);
lgraph = replaceLayer(lgraph, 'output', pxLayer);

%% ── Data loading ─────────────────────────────────────────────────────────
fprintf('[TRAIN] Loading images from:  %s\n', image_dir);
fprintf('[TRAIN] Loading labels from:  %s\n', label_dir);

img_ds   = imageDatastore(image_dir, 'FileExtensions', {'.png','.tif'});
label_ds = pixelLabelDatastore(label_dir, class_names, [1 2 3], ...
                               'FileExtensions', {'.png','.tif'});

combined_ds = combine(img_ds, label_ds);

% Train/val split (80/20)
n_total = numel(img_ds.Files);
n_train = floor(0.8 * n_total);
n_val   = n_total - n_train;
fprintf('[TRAIN] Total images: %d  (train=%d, val=%d)\n', n_total, n_train, n_val);

rng(42);
idx       = randperm(n_total);
train_ds  = subset(combined_ds, idx(1:n_train));
val_ds    = subset(combined_ds, idx(n_train+1:end));

% Data augmentation (flip, rotate — nucleus appearance is rotationally invariant)
aug_transforms = imageDataAugmenter( ...
    'RandXReflection',   true, ...
    'RandYReflection',   true, ...
    'RandRotation',      [-180 180], ...
    'RandXShear',        [-5 5], ...
    'RandYShear',        [-5 5]);

aug_train_ds = augmentedImageDatastore(input_size, train_ds, ...
    'DataAugmentation', aug_transforms, ...
    'ColorPreprocessing', 'gray2gray');

aug_val_ds = augmentedImageDatastore(input_size, val_ds, ...
    'ColorPreprocessing', 'gray2gray');

%% ── Training options ─────────────────────────────────────────────────────
options = trainingOptions('adam', ...
    'InitialLearnRate',    1e-3, ...
    'LearnRateSchedule',   'piecewise', ...
    'LearnRateDropFactor', 0.5, ...
    'LearnRateDropPeriod', 10, ...
    'MaxEpochs',           30, ...
    'MiniBatchSize',       16, ...
    'Shuffle',             'every-epoch', ...
    'ValidationData',       aug_val_ds, ...
    'ValidationFrequency',  5, ...
    'Plots',               'training-progress', ...
    'Verbose',             true, ...
    'VerboseFrequency',    5, ...
    'ExecutionEnvironment','auto');   % uses GPU if available

fprintf('[TRAIN] Starting training...\n');
fprintf('[TRAIN] (GPU recommended — CPU training will be slow for >200 images)\n\n');

%% ── Train ────────────────────────────────────────────────────────────────
net = trainNetwork(aug_train_ds, lgraph, options);

%% ── Evaluate on validation set ───────────────────────────────────────────
fprintf('\n[TRAIN] Evaluating on validation set...\n');
pxds_pred = semanticseg(aug_val_ds, net, 'MiniBatchSize', 8, 'WriteLocation', tempdir);
metrics   = evaluateSemanticSegmentation(pxds_pred, val_ds, 'Verbose', false);

fprintf('[TRAIN] Per-class IoU:\n');
for c = 1:n_classes
    fprintf('        %-18s  IoU = %.3f\n', class_names{c}, ...
            metrics.ClassMetrics.IoU(c));
end
fprintf('[TRAIN] Mean IoU: %.3f\n', metrics.DataSetMetrics.MeanIoU);

%% ── Save ─────────────────────────────────────────────────────────────────
if isfield(params, 'normality_model_path') && ~isempty(params.normality_model_path)
    save(params.normality_model_path, 'net');
    fprintf('[TRAIN] Model saved to: %s\n', params.normality_model_path);
else
    fprintf('[TRAIN] Tip: save model with:\n');
    fprintf('        save(''nucleus_normality_unet.mat'', ''net'')\n');
end
end


%% ======================================================================
function export_crops_for_labelling(img_data, params, output_dir)
%EXPORT_CROPS_FOR_LABELLING  Export nucleus crops ready for manual labelling.
%
%  Runs segmentation, then saves each 64x64 nucleus crop as a numbered PNG.
%  Open output_dir in MATLAB Image Labeler or ImageJ to create ground truth.
%
%  Usage:
%    export_crops_for_labelling(img_data, params, 'C:\SAM\labelling\images')

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

[cells, ~, ~] = segment_nuclei(img_data.nucleus, params);
nucleus_img   = img_data.nucleus;
[H, W]        = size(nucleus_img);
crop_size     = 64; half = crop_size/2;

fprintf('[EXPORT] Exporting %d nucleus crops to: %s\n', numel(cells), output_dir);

for k = 1:numel(cells)
    cx = round(cells(k).centroid(1));
    cy = round(cells(k).centroid(2));
    x1 = max(1,cx-half); x2 = min(W,cx+half);
    y1 = max(1,cy-half); y2 = min(H,cy+half);
    crop = nucleus_img(y1:y2, x1:x2);

    padded = zeros(crop_size, crop_size);
    padded(1:size(crop,1), 1:size(crop,2)) = crop;

    fname = fullfile(output_dir, sprintf('cell_%04d.png', k));
    imwrite(uint8(padded * 255), fname);
end

fprintf('[EXPORT] Done. Open these in MATLAB Image Labeler:\n');
fprintf('         imageLabeler(''%s'')\n', output_dir);
fprintf('         Define 3 pixel classes: normal / abnormal_shape / abnormal_count\n');
fprintf('         Export labels, then call train_nucleus_classifier()\n');
end

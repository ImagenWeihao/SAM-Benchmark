function [cells_out, model_info] = unet_cellpose(cells_in, img_data, params)
%UNET_CELLPOSE  Cellpose adapter using pretrained 'nuclei' model.
%
%  Uses Cellpose 2.0 pretrained 'nuclei' model — no training needed.
%  Derives normality score from Cellpose flow fields and mask properties:
%    - Flow coherence  : how well the predicted flow converges (shape quality)
%    - Mask coverage   : fraction of crop covered by detected nucleus
%    - n_masks         : number of separate objects detected in crop
%
%  Requirements:
%    pip install cellpose
%    pyenv('Version', 'C:\path\to\python.exe')  % in MATLAB once per session
%
%  params fields used:
%    .cellpose_model   — model type: 'nuclei' (default) | 'cyto' | 'cyto2'
%    .nucleus_channel  — channel index for DAPI image
%    .unet_threshold   — flow threshold for Cellpose (default: 0.4)

%% ── Model info ───────────────────────────────────────────────────────────
cp_model_type = 'nuclei';
if isfield(params, 'cellpose_model') && ~isempty(params.cellpose_model)
    cp_model_type = params.cellpose_model;
end

model_info.name    = 'Cellpose';
model_info.version = cp_model_type;
model_info.target  = params.unet_target;
model_info.status  = 'pending';

cells_out = cells_in;
n         = numel(cells_in);

%% ── Check Python + Cellpose ──────────────────────────────────────────────
try
    py.importlib.import_module('cellpose');
catch ME
    error('UNET:cellposeNotFound', ...
          ['Cellpose not found in Python environment.\n' ...
           'Fix:\n' ...
           '  1. pip install cellpose\n' ...
           '  2. In MATLAB: pyenv(''Version'', ''C:\\path\\to\\python.exe'')\n' ...
           'Original error: %s'], ME.message);
end

%% ── Load model ───────────────────────────────────────────────────────────
fprintf('[Cellpose] Loading pretrained model: %s\n', cp_model_type);
cp_models = py.importlib.import_module('cellpose.models');
np        = py.importlib.import_module('numpy');
model     = cp_models.Cellpose(pyargs('model_type', cp_model_type, 'gpu', false));

%% ── Population-level area stats ─────────────────────────────────────────
all_areas = [cells_in.area];
mean_area = mean(all_areas);
std_area  = std(all_areas);

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

    % Scale to uint16 range — Cellpose handles this well for DAPI
    crop_u16 = uint16(crop ./ max(crop(:)+eps) * 65535);
    py_img   = np.array(crop_u16);
    py_imgs  = py.list({py_img});

    try
        % channels=[0,0] = grayscale input
        result = model.eval(py_imgs, ...
                     pyargs('diameter',   py.None, ...
                            'channels',   py.list({py.int(0), py.int(0)}), ...
                            'flow_threshold', params.unet_threshold));

        masks   = result{1};   % list of label arrays
        flows   = result{2};   % list of flow arrays [dy, dx, cellprob, dist]
        styles  = result{3};   % style vector (optional)

        mask_mat  = double(masks{1});
        n_masks   = double(max(mask_mat(:)));

        % Flow coherence: use cell probability map (flows{1}{2} = cellprob)
        flow_arr  = result{2}{1};   % flows for first image
        cellprob  = double(py.numpy.array(flow_arr{3}));  % index 3 = cellprob
        mean_prob = mean(cellprob(:));

        % Mask coverage: fraction of crop with detected nucleus
        coverage  = sum(mask_mat(:) > 0) / numel(mask_mat);

    catch ME
        warning('[Cellpose] Cell %d inference failed: %s', k, ME.message);
        cells_out(k).unet_score  = 0.5;
        cells_out(k).unet_class  = 'normal';
        cells_out(k).cp_n_masks  = 0;
        cells_out(k).cp_prob     = 0;
        cells_out(k).cp_coverage = 0;
        continue;
    end

    %% Size anomaly
    area_z = abs(cells_in(k).area - mean_area) / (std_area + eps);

    %% Assign class
    if n_masks >= 2
        % Multiple masks → binucleated or micronuclei (C, E)
        unet_class = 'abnormal_count';
    elseif area_z > 2.0
        % Size outlier → polyploid or micronucleus (D, E)
        unet_class = 'abnormal_count';
    elseif coverage < 0.15 || mean_prob < 0.3
        % Poor detection → likely irregular shape (B, F)
        unet_class = 'abnormal_shape';
    else
        unet_class = 'normal';
    end

    %% Score = combination of prob, coverage, size normality
    size_penalty = min(1.0, area_z * 0.3);
    score = mean_prob * coverage * (1 - size_penalty);
    score = max(0, min(1, score));

    cells_out(k).unet_score  = score;
    cells_out(k).unet_class  = unet_class;
    cells_out(k).cp_n_masks  = n_masks;
    cells_out(k).cp_prob     = mean_prob;
    cells_out(k).cp_coverage = coverage;

    fprintf('[Cellpose] Cell %2d: n_masks=%d  prob=%.3f  coverage=%.3f  → %s (score=%.3f)\n', ...
            cells_in(k).id, n_masks, mean_prob, coverage, unet_class, score);
end

model_info.status      = 'ok';
model_info.model_name  = cp_model_type;
model_info.n_processed = n;
end

function [cells_out, model_info] = unet_matlab(cells_in, img_data, params)
%UNET_MATLAB  MATLAB Deep Learning Toolbox UNet adapter.
%
%  Requires:
%    - Deep Learning Toolbox
%    - A trained .mat file containing a 'net' variable
%      Set params.unet_model_path to its location.
%
%  If no model file is provided, falls back to mock scoring with a warning.
%  To train your own model use train_nucleus_classifier.m.

model_info.name   = 'MATLAB UNet';
model_info.target = params.unet_target;
cells_out         = cells_in;
n                 = numel(cells_in);

% Try to load trained model
model_loaded = false;
if isfield(params, 'unet_model_path') && ~isempty(params.unet_model_path) ...
        && isfile(params.unet_model_path)
    fprintf('[UNET] Loading MATLAB UNet from: %s\n', params.unet_model_path);
    loaded       = load(params.unet_model_path, 'net');
    net          = loaded.net;
    model_loaded = true;
    model_info.version = 'custom trained';
else
    warning('[UNET] No model file at params.unet_model_path. Falling back to mock scoring.');
    fprintf('[UNET] To train: call train_nucleus_classifier() and set params.unet_model_path.\n');
    model_info.version = 'untrained — mock fallback';
end

nucleus_img = img_data.nucleus;
[H, W]      = size(nucleus_img);
crop_size   = 128;
half        = crop_size / 2;

if model_loaded
    for k = 1:n
        cx = round(cells_in(k).centroid(1));
        cy = round(cells_in(k).centroid(2));
        x1 = max(1, cx-half); x2 = min(W, cx+half);
        y1 = max(1, cy-half); y2 = min(H, cy+half);
        crop = nucleus_img(y1:y2, x1:x2);

        % Pad to crop_size x crop_size
        padded = zeros(crop_size, crop_size);
        padded(1:size(crop,1), 1:size(crop,2)) = crop;

        try
            seg_map = semanticseg(single(padded), net);
            counts  = histcounts(uint8(seg_map), 1:4);
            [conf, cls] = max(counts);
            cells_out(k).unet_score = double(conf) / sum(counts);
            cells_out(k).unet_class = params.unet_target;
        catch
            cells_out(k).unet_score = 0.5;
            cells_out(k).unet_class = params.unet_target;
        end
    end
else
    % Fall back to mock
    [cells_out, ~] = unet_mock(cells_in, img_data, params);
    for k = 1:n
        cells_out(k).unet_class = ['[matlab-fallback] ' cells_out(k).unet_class];
    end
end

model_info.status      = 'ok';
model_info.n_processed = n;
end

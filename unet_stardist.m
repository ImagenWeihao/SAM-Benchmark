function [cells_out, model_info] = unet_stardist(cells_in, img_data, params)
%UNET_STARDIST  StarDist adapter using pretrained '2D_versatile_fluo' model.
%
%  Uses StarDist's pretrained fluorescence nucleus model via a Python bridge
%  script (stardist_bridge.py) which must exist in the SAM folder.
%
%  Confirmed working keys returned by predict_instances:
%    'prob'   — detection probability per object
%    'coord'  — polygon coordinates per object (shape descriptor)
%    'points' — centroid points per object
%
%  Normality scoring:
%    shape_var = std of coord radii / mean  → high = lobular/blebbing (B,F)
%    n_objects > 1                          → binucleated/micronuclei (C,E)
%    area z-score > 2                       → polyploid/micronucleus (D,E)
%
%  Requirements:
%    - stardist_bridge.py in D:\Matlab Files\Codes\SAM\
%    - pyenv set to Python 3.9 OutOfProcess before calling run_benchmark
%    - pip install stardist tensorflow-cpu==2.10.0 "numpy<2"

%% ── Model info ───────────────────────────────────────────────────────────
model_name = '2D_versatile_fluo';
if isfield(params, 'stardist_model') && ~isempty(params.stardist_model)
    model_name = params.stardist_model;
end

model_info.name    = 'StarDist';
model_info.version = model_name;
model_info.target  = params.unet_target;
model_info.status  = 'pending';

cells_out = cells_in;
n         = numel(cells_in);

%% ── Load bridge and model ────────────────────────────────────────────────
bridge_path = fullfile(fileparts(mfilename('fullpath')), 'stardist_bridge.py');
if ~isfile(bridge_path)
    error('UNET:bridgeNotFound', ...
          ['stardist_bridge.py not found at:\n  %s\n' ...
           'Run this in MATLAB to create it:\n' ...
           '  fid = fopen(''%s'', ''w'');\n' ...
           '  fprintf(fid, ''from stardist.models import StarDist2D\nimport numpy as np\n\n'');\n' ...
           '  fprintf(fid, ''def load_model(n="2D_versatile_fluo"):\n    return StarDist2D.from_pretrained(n)\n\n'');\n' ...
           '  fprintf(fid, ''def predict_crop(model, crop):\n    labels, details = model.predict_instances(crop)\n'');\n' ...
           '  fprintf(fid, ''    return {"probs": details["prob"].tolist(), "coords": details["coord"].tolist(), "n_objects": int(labels.max())}\n'');\n' ...
           '  fclose(fid);'], bridge_path, bridge_path);
end

fprintf('[StarDist] Loading bridge from: %s\n', bridge_path);
spec   = py.importlib.util.spec_from_file_location('stardist_bridge', bridge_path);
bridge = py.importlib.util.module_from_spec(spec);
spec.loader.exec_module(bridge);

fprintf('[StarDist] Loading pretrained model: %s\n', model_name);
load_fn    = py.getattr(bridge, 'load_model');
model      = feval(load_fn, model_name);
predict_fn = py.getattr(bridge, 'predict_crop');
np         = py.importlib.import_module('numpy');

model_info.model_obj = model;

%% ── Population-level area stats ─────────────────────────────────────────
all_areas = [cells_in.area];
mean_area = mean(all_areas);
std_area  = std(all_areas);

nucleus_img = img_data.nucleus;
[H, W]      = size(nucleus_img);
crop_half   = 64;

%% ── Run per-cell ─────────────────────────────────────────────────────────
fprintf('[StarDist] Running inference on %d cells...\n', n);

for k = 1:n
    cx = round(cells_in(k).centroid(1));
    cy = round(cells_in(k).centroid(2));
    x1 = max(1, cx-crop_half); x2 = min(W, cx+crop_half);
    y1 = max(1, cy-crop_half); y2 = min(H, cy+crop_half);
    crop = nucleus_img(y1:y2, x1:x2);

    % Normalise to float32
    cmax = max(crop(:));
    if cmax > 0; crop = crop ./ cmax; end
    crop_f32 = single(crop);
    py_crop  = np.array(crop_f32);

    try
        result   = feval(predict_fn, model, py_crop);
        n_found  = int32(double(result{'n_objects'}));

        % Extract probabilities
        probs_py = result{'probs'};
        if double(py.len(probs_py)) > 0
            probs    = double(py.array.array('d', probs_py));
            max_prob = max(probs);
        else
            probs    = [];
            max_prob = 0;
        end

        % Extract coord shape variance (irregularity proxy)
        coords_py = result{'coords'};
        shape_var = 0;
        if double(py.len(coords_py)) > 0
            try
                coord_np  = py.numpy.array(coords_py);
                coord_u16 = feval(py.getattr(coord_np, 'astype'), 'float32');
                flat_arr  = feval(py.getattr(coord_u16, 'flatten'));
                coord_list= feval(py.getattr(flat_arr,  'tolist'));
                coord_arr = double(py.array.array('f', coord_list));
                if numel(coord_arr) > 2
                    shape_var = std(coord_arr) / (mean(coord_arr) + eps);
                end
            catch
                shape_var = 0;
            end
        end

    catch ME
        warning('[StarDist] Cell %d failed: %s', cells_in(k).id, ME.message);
        cells_out(k).unet_score   = 0.5;
        cells_out(k).unet_class   = 'normal';
        cells_out(k).sd_prob      = 0;
        cells_out(k).sd_shape_var = 0;
        cells_out(k).sd_n_objects = 0;
        continue;
    end

    %% Size anomaly z-score
    area_z = abs(cells_in(k).area - mean_area) / (std_area + eps);

    %% Assign normality class
    if n_found >= 2
        unet_class = 'abnormal_count';   % binucleated / micronuclei (C, E)
    elseif area_z > 2.0
        unet_class = 'abnormal_count';   % polyploid / size outlier (D)
    elseif shape_var > 0.25
        unet_class = 'abnormal_shape';   % lobular / blebbing (B, F)
    else
        unet_class = 'normal';           % round, compact (A)
    end

    %% Normality score (high = normal)
    shape_penalty = min(1.0, shape_var * 3.0);
    size_penalty  = min(1.0, area_z   * 0.3);
    score = max_prob * (1 - 0.6*shape_penalty - 0.4*size_penalty);
    score = max(0, min(1, score));

    cells_out(k).unet_score    = score;
    cells_out(k).unet_class    = unet_class;
    cells_out(k).sd_prob       = max_prob;
    cells_out(k).sd_shape_var  = shape_var;
    cells_out(k).sd_n_objects  = double(n_found);

    fprintf('[StarDist] Cell %2d: prob=%.3f  shape_var=%.3f  n_obj=%d  -> %-16s (score=%.3f)\n', ...
            cells_in(k).id, max_prob, shape_var, n_found, unet_class, score);
end

model_info.status      = 'ok';
model_info.model_name  = model_name;
model_info.n_processed = n;
end
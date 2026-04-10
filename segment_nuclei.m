function [cells, seg_mask, bw_nuclei] = segment_nuclei(nucleus_img, params)
%SEGMENT_NUCLEI  Detect and measure all nuclei in a nucleus-channel image.
%
%  [cells, seg_mask, bw_nuclei] = segment_nuclei(nucleus_img, params)
%
%  Inputs:
%    nucleus_img  — normalised [0,1] double grayscale image (DAPI / 405 nm)
%    params       — benchmark params struct
%
%  Outputs:
%    cells        — struct array, one entry per detected nucleus:
%                     .id           unique integer ID
%                     .centroid     [x, y] in pixels
%                     .area         area in pixels²
%                     .circularity  4π·Area / Perimeter²  ∈ (0,1]
%                     .bbox         [x y w h] bounding box
%                     .mean_int     mean nucleus intensity
%                     .perimeter    perimeter in pixels
%    seg_mask     — label matrix (uint16), 0=background, N=cell N
%    bw_nuclei    — binary mask before labelling

fprintf('[SEG]  Segmenting nuclei...\n');

%% ── Pre-processing ───────────────────────────────────────────────────────
% Mild Gaussian smoothing to reduce shot noise
img_smooth = imgaussfilt(nucleus_img, 1.5);

% Adaptive thresholding (handles uneven illumination across well)
sensitivity = 0.5;
bw_adapt    = imbinarize(img_smooth, 'adaptive', ...
                         'Sensitivity', sensitivity, ...
                         'ForegroundPolarity', 'bright');

% Morphological cleanup
se_open  = strel('disk', 3);
se_close = strel('disk', 5);
bw_clean = imopen(bw_adapt,  se_open);
bw_clean = imclose(bw_clean, se_close);
bw_clean = imfill(bw_clean, 'holes');

% Remove very small objects (debris) and very large (clusters / artefacts)
min_area_px = 200;    % ~15 µm diameter at 0.3 µm/px
max_area_px = 80000;  % filter out large aggregates
bw_clean    = bwareaopen(bw_clean, min_area_px);
bw_clean    = bw_clean & ~bwareaopen(bw_clean, max_area_px);

% Separate touching nuclei with watershed
dist_map  = -bwdist(~bw_clean);
dist_map  = imhmin(dist_map, 2);   % suppress shallow minima
ws_labels = watershed(dist_map);
bw_nuclei = bw_clean & (ws_labels ~= 0);

%% ── Label and measure ────────────────────────────────────────────────────
seg_mask = uint16(bwlabel(bw_nuclei));
props    = regionprops(seg_mask, nucleus_img, ...
    'Area', 'Centroid', 'Perimeter', 'BoundingBox', 'MeanIntensity');

n_detected = numel(props);
fprintf('[SEG]  Detected %d nuclei.\n', n_detected);

if n_detected == 0
    cells = struct([]);
    return;
end

%% ── Compute circularity and package struct array ─────────────────────────
cells = struct( ...
    'id',          num2cell((1:n_detected)'), ...
    'centroid',    {props.Centroid}', ...
    'area',        num2cell([props.Area]'), ...
    'perimeter',   num2cell([props.Perimeter]'), ...
    'circularity', num2cell(zeros(n_detected,1)), ...
    'bbox',        {props.BoundingBox}', ...
    'mean_int',    num2cell([props.MeanIntensity]'));

for i = 1:n_detected
    p = cells(i).perimeter;
    a = cells(i).area;
    if p > 0
        cells(i).circularity = (4 * pi * a) / (p^2);
    else
        cells(i).circularity = 0;
    end
end

% Clamp circularity to [0,1] (floating point can give tiny overshoots)
circ_vals = [cells.circularity];
circ_vals = min(circ_vals, 1.0);
for i = 1:n_detected
    cells(i).circularity = circ_vals(i);
end

fprintf('[SEG]  Circularity range: [%.3f – %.3f]  mean=%.3f\n', ...
        min(circ_vals), max(circ_vals), mean(circ_vals));
end

function [cells, seg_mask, bw_nuclei] = segment_nuclei_l1(nucleus_img, params) %#ok<INUSD>
%SEGMENT_NUCLEI_L1  Otsu + connected-components baseline (Level 1).
%
%  Intentionally minimal: no morphological cleanup, no area filter, no
%  watershed splitting, no circularity rule. Used as the lower-tier
%  benchmark against L2 (rule-based pipeline + circularity filter) and
%  L3 (StarDist ML model).
%
%  [cells, seg_mask, bw_nuclei] = segment_nuclei_l1(nucleus_img, params)
%
%  Returns the same struct shape as segment_nuclei.m so downstream code is
%  drop-in compatible. params is accepted for signature parity but unused.
fprintf('[SEG-L1] Otsu + connected components...\n');
%% -- Global Otsu threshold -----------------------------------------------
% graythresh picks the threshold that minimises intra-class variance over
% the full intensity histogram. Data-driven, no if-then logic.
level     = graythresh(nucleus_img);
bw_nuclei = imbinarize(nucleus_img, level);
%% -- Label connected components ------------------------------------------
seg_mask = uint16(bwlabel(bw_nuclei));
props    = regionprops(seg_mask, nucleus_img, ...
    'Area', 'Centroid', 'Perimeter', 'BoundingBox', 'MeanIntensity', 'Solidity');
n_detected = numel(props);
fprintf('[SEG-L1] Detected %d connected components.\n', n_detected);
if n_detected == 0
    cells = struct([]);
    return;
end
%% -- Package cells struct (matches segment_nuclei output schema) ---------
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

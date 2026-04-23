function [cells_out, model_info] = classify_cells_unet(cells_in, img_data, params)
%CLASSIFY_CELLS_UNET  Dispatch cell crops to the selected UNet classifier.
%
%  [cells_out, model_info] = classify_cells_unet(cells_in, img_data, params)
%
%  Adds .unet_score (0-1) and .unet_class fields to each cell.
%  Dispatches to standalone adapter files based on params.unet_model:
%
%    'stardist'    -- unet_stardist.m    pip install stardist tensorflow-cpu
%    'cellpose'    -- unet_cellpose.m    pip install cellpose
%    'nuclear_seg' -- unet_nuclear_seg.m pip install nuclear-segmentator
%    'matlab'      -- unet_matlab.m      Deep Learning Toolbox + trained .mat
%    'mock'        -- unet_mock.m        No install needed (rule-based proxy)
fprintf('[UNET] Model: %s  |  Target: "%s"\n', params.unet_model, params.unet_target);
fprintf('[UNET] Classifying %d cells...\n', numel(cells_in));
t0 = tic;
switch lower(params.unet_model)
    case 'stardist'
        [cells_out, model_info] = unet_stardist(cells_in, img_data, params);
    case 'cellpose'
        [cells_out, model_info] = unet_cellpose(cells_in, img_data, params);
    case 'nuclear_seg'
        [cells_out, model_info] = unet_nuclear_seg(cells_in, img_data, params);
    case 'matlab'
        [cells_out, model_info] = unet_matlab(cells_in, img_data, params);
    case 'mock'
        [cells_out, model_info] = unet_mock(cells_in, img_data, params);
    otherwise
        error('UNET:unknownModel', ...
              ['Unknown model: "%s"\n' ...
               'Valid options:\n' ...
               '  ''stardist''    -- pip install stardist tensorflow-cpu==2.10.0\n' ...
               '  ''cellpose''    -- pip install cellpose\n' ...
               '  ''nuclear_seg'' -- pip install nuclear-segmentator\n' ...
               '  ''matlab''      -- Deep Learning Toolbox + trained .mat\n' ...
               '  ''mock''        -- no install needed\n'], params.unet_model);
end
elapsed = toc(t0);
scores  = [cells_out.unet_score];
fprintf('[UNET] Done in %.2f s | min=%.3f  max=%.3f  mean=%.3f\n', ...
        elapsed, min(scores), max(scores), mean(scores));
end

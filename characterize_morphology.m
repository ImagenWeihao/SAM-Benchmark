function morph = characterize_morphology(target_cells, zstack_results, img_data, params)
%CHARACTERIZE_MORPHOLOGY  Cell morphology characterization (placeholder).
%
%  morph = characterize_morphology(target_cells, zstack_results, img_data, params)
%
%  Currently computes a standard set of measurable features from the
%  segmented nuclei and Z-stack data. This function is designed as a
%  placeholder -- additional characterization pipelines (e.g., 3D shape
%  reconstruction, multi-channel co-localisation, texture features) can be
%  added here in future levels.
%
%  Output struct array per cell:
%    .cell_id          -- cell index
%    .nucleus_area     -- area in pixels^2
%    .circularity      -- 4??A/P^2
%    .mean_intensity   -- mean nucleus channel intensity
%    .texture_contrast -- GLCM contrast (grayscale co-occurrence texture)
%    .aspect_ratio     -- major/minor axis length ratio (elongation)
%    .solidity         -- area / convex hull area
%    .focus_quality    -- best-focus Laplacian variance score
fprintf('[MORPH] Characterising morphology for %d cells...\n', numel(target_cells));
nucleus_img = img_data.nucleus;
n_cells     = numel(target_cells);
morph       = struct( ...
    'cell_id',          num2cell(zeros(n_cells,1)), ...
    'nucleus_area',     num2cell(zeros(n_cells,1)), ...
    'circularity',      num2cell(zeros(n_cells,1)), ...
    'mean_intensity',   num2cell(zeros(n_cells,1)), ...
    'texture_contrast', num2cell(zeros(n_cells,1)), ...
    'aspect_ratio',     num2cell(zeros(n_cells,1)), ...
    'solidity',         num2cell(zeros(n_cells,1)), ...
    'focus_quality',    num2cell(zeros(n_cells,1)));
for k = 1:n_cells
    cell_k = target_cells(k);
    %% Basic nucleus features (from segmentation)
    morph(k).cell_id        = cell_k.id;
    morph(k).nucleus_area   = cell_k.area;
    morph(k).circularity    = cell_k.circularity;
    morph(k).mean_intensity = cell_k.mean_int;
    %% Aspect ratio and solidity from regionprops on nucleus crop
    bbox = round(cell_k.bbox);     % [x y w h]
    x1 = max(1, bbox(1));
    y1 = max(1, bbox(2));
    x2 = min(size(nucleus_img,2), bbox(1)+bbox(3));
    y2 = min(size(nucleus_img,1), bbox(2)+bbox(4));
    crop = nucleus_img(y1:y2, x1:x2);
    bw_crop = imbinarize(crop, 'adaptive', 'Sensitivity', 0.5);
    rp      = regionprops(bw_crop, 'MajorAxisLength', 'MinorAxisLength', 'Solidity');
    if ~isempty(rp)
        [~, largest] = max([rp.MajorAxisLength]);
        if rp(largest).MinorAxisLength > 0
            morph(k).aspect_ratio = rp(largest).MajorAxisLength / rp(largest).MinorAxisLength;
        else
            morph(k).aspect_ratio = 1;
        end
        morph(k).solidity = rp(largest).Solidity;
    else
        morph(k).aspect_ratio = 1;
        morph(k).solidity     = 1;
    end
    %% Texture: GLCM contrast on best-focus Z plane
    zr = zstack_results(k);
    if ~isempty(zr.z_planes)
        best_plane  = zr.z_planes{zr.best_focus_z};
        img_uint8   = uint8(best_plane * 255);
        glcm        = graycomatrix(img_uint8, 'Offset', [0 1], 'Symmetric', true);
        stats       = graycoprops(glcm, 'Contrast');
        morph(k).texture_contrast = stats.Contrast;
        morph(k).focus_quality    = max(zr.focus_scores);
    end
end
fprintf('[MORPH] Characterisation complete.\n');
%% Print summary table
fprintf('\n%-6s  %-10s  %-12s  %-13s  %-13s  %-9s\n', ...
        'CellID','Area(px^2)','Circularity','AspectRatio','Solidity','FocusQual');
fprintf('%s\n', repmat('-',1,70));
for k = 1:n_cells
    fprintf('%-6d  %-10.0f  %-12.3f  %-13.3f  %-13.3f  %-9.4f\n', ...
            morph(k).cell_id, morph(k).nucleus_area, morph(k).circularity, ...
            morph(k).aspect_ratio, morph(k).solidity, morph(k).focus_quality);
end
fprintf('\n');
end

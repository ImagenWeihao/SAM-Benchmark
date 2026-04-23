function zstack_results = acquire_zstack(target_cells, img_data, params)
%ACQUIRE_ZSTACK  Simulate or extract Z-stack acquisition for selected cells.
%
%  zstack_results = acquire_zstack(target_cells, img_data, params)
%
%  If the .nd2 file contains real Z planes (img_data.n_z > 1), this function
%  extracts the actual Z-stack for each target cell's bounding box.
%  If n_z == 1, it simulates a Z-stack by applying synthetic defocus blur
%  to mimic out-of-focus planes above and below the in-focus plane.
%
%  Output struct array per cell:
%    .cell_id      -- original cell ID from segmentation
%    .centroid     -- [x, y]
%    .z_planes     -- {n_z +/- 1} cell of cropped images per Z level
%    .best_focus_z -- estimated in-focus plane index (Laplacian variance)
%    .focus_scores -- focus quality score per Z plane
n_targets = numel(target_cells);
fprintf('[ZSTK] Acquiring Z-stacks for %d target cells (%d Z-planes each)...\n', ...
        n_targets, params.z_planes);
nucleus_img = img_data.nucleus;
[H, W]      = size(nucleus_img);
% Half-size of crop window around each nucleus centroid (pixels)
crop_half = 64;
zstack_results = struct( ...
    'cell_id',      num2cell(zeros(n_targets,1)), ...
    'centroid',     cell(n_targets,1), ...
    'z_planes',     cell(n_targets,1), ...
    'best_focus_z', num2cell(zeros(n_targets,1)), ...
    'focus_scores', cell(n_targets,1));
for k = 1:n_targets
    cell_k = target_cells(k);
    cx = round(cell_k.centroid(1));
    cy = round(cell_k.centroid(2));
    % Clamp crop to image boundaries
    x1 = max(1, cx - crop_half);
    x2 = min(W, cx + crop_half);
    y1 = max(1, cy - crop_half);
    y2 = min(H, cy + crop_half);
    zstack_results(k).cell_id  = cell_k.id;
    zstack_results(k).centroid = cell_k.centroid;
    if img_data.n_z > 1
        %% Real Z-stack extraction
        n_planes = min(params.z_planes, img_data.n_z);
        z_imgs   = cell(n_planes, 1);
        raw_stack = img_data.raw{1,1};  % full plane list from bfopen
        n_ch      = img_data.n_channels;
        for z = 1:n_planes
            % Plane index: first T, nucleus channel c, plane z
            % bfopen order: C varies fastest -> plane = (z-1)*n_ch + c
            plane_idx = (z-1)*n_ch + params.nucleus_channel;
            if plane_idx <= size(raw_stack, 1)
                raw_plane = double(raw_stack{plane_idx, 1});
                cmax = max(raw_plane(:));
                if cmax > 0; raw_plane = raw_plane ./ cmax; end
                z_imgs{z} = raw_plane(y1:y2, x1:x2);
            else
                z_imgs{z} = nucleus_img(y1:y2, x1:x2);
            end
        end
    else
        %% Simulated Z-stack using synthetic defocus
        n_planes   = params.z_planes;
        z_imgs     = cell(n_planes, 1);
        in_focus   = nucleus_img(y1:y2, x1:x2);
        mid        = ceil(n_planes / 2);
        for z = 1:n_planes
            defocus_steps = abs(z - mid);         % 0 = in-focus plane
            sigma = defocus_steps * 1.2;           % blur increases with distance
            if sigma > 0
                z_imgs{z} = imgaussfilt(in_focus, sigma);
            else
                z_imgs{z} = in_focus;
            end
            % Add realistic photon noise that increases off-focus
            noise_level = 0.01 + defocus_steps * 0.005;
            z_imgs{z}   = z_imgs{z} + noise_level * randn(size(z_imgs{z}));
            z_imgs{z}   = max(0, min(1, z_imgs{z}));
        end
    end
    %% Focus quality: Laplacian variance (higher = sharper)
    focus_scores = zeros(1, n_planes);
    for z = 1:n_planes
        lap = del2(z_imgs{z});
        focus_scores(z) = var(lap(:));
    end
    [~, best_z] = max(focus_scores);
    zstack_results(k).z_planes     = z_imgs;
    zstack_results(k).best_focus_z = best_z;
    zstack_results(k).focus_scores = focus_scores;
end
fprintf('[ZSTK] Z-stack acquisition complete.\n');
end

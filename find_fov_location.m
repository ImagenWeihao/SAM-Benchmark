function result = find_fov_location(path_10x, path_40x, params)
%FIND_FOV_LOCATION  Find the location of a 40X FOV within a 10X overview image
%                   using normalised cross-correlation.
%
%  Workflow:
%    1. Load 10X image (single Z, 405nm channel)
%    2. Load 40X image and compute MIP across all Z planes (Series 1-15)
%    3. Estimate physical scale factor from pixel sizes
%    4. Resize 40X MIP to match 10X pixel scale
%    5. Run normalised cross-correlation (normxcorr2)
%    6. Find peak location -> convert to physical coordinates
%    7. Display result
%
%  Usage:
%    params.channel_10x   = 1;     % 405nm channel index in 10X file
%    params.channel_40x   = 2;     % 405nm channel index in 40X file
%    params.series_40x    = 1;     % first series of FOV1 in 40X file
%    params.series_40x_end= 15;    % last series of FOV1 in 40X file
%    result = find_fov_location(path_10x, path_40x, params);
if nargin < 3; params = struct(); end
if ~isfield(params, 'channel_10x');    params.channel_10x    = 1;  end
if ~isfield(params, 'channel_40x');    params.channel_40x    = 2;  end
if ~isfield(params, 'series_40x');     params.series_40x     = 1;  end
if ~isfield(params, 'series_40x_end'); params.series_40x_end = 15; end
fprintf('\n================================================\n');
fprintf('  FOV LOCATION FINDER -- 40X in 10X (all FOVs)\n');
fprintf('================================================\n\n');
%% -- Step 1: Load ALL 10X series -----------------------------------------
fprintf('[10X]  Loading overview image...\n');
fprintf('       %s\n', path_10x);
data_10x = bfopen(path_10x);
n_series_10x = size(data_10x, 1);
fprintf('[10X]  Loaded %d series (FOVs).\n', n_series_10x);
omeMeta_10x = data_10x{1,4};
try
    px_10x = double(omeMeta_10x.getPixelsPhysicalSizeX(0).value());
catch
    px_10x = 0.65;
    warning('[10X]  Could not read pixel size, assuming %.4f um/px', px_10x);
end
fprintf('[10X]  Pixel size: %.4f um/px\n', px_10x);
% Load 405nm channel from each series
n_planes_10x = size(data_10x{1,1}, 1);
fprintf('[10X]  Planes per series: %d\n', n_planes_10x);
imgs_10x = cell(1, n_series_10x);
for s = 1:n_series_10x
    raw = double(data_10x{s,1}{params.channel_10x, 1});
    imgs_10x{s} = raw ./ (max(raw(:)) + eps);
end
[H10, W10] = size(imgs_10x{1});
fprintf('[10X]  Image size per FOV: %d x %d px  (%.1f x %.1f um)\n', ...
        W10, H10, W10*px_10x, H10*px_10x);
%% -- Step 2: Load 40X MIP -------------------------------------------------
fprintf('\n[40X]  Loading and computing MIP (Series %d-%d)...\n', ...
        params.series_40x, params.series_40x_end);
fprintf('       %s\n', path_40x);
data_40x = bfopen(path_40x);
fprintf('[40X]  Loaded %d series total.\n', size(data_40x,1));
omeMeta_40x = data_40x{1,4};
try
    px_40x = double(omeMeta_40x.getPixelsPhysicalSizeX(0).value());
catch
    px_40x = 0.1625;
    warning('[40X]  Could not read pixel size, assuming %.4f um/px', px_40x);
end
fprintf('[40X]  Pixel size: %.4f um/px\n', px_40x);
% Compute MIP across series range
img_40x_mip = zeros(size(data_40x{params.series_40x,1}{1,1}));
n_z = params.series_40x_end - params.series_40x + 1;
for s = params.series_40x:params.series_40x_end
    plane = double(data_40x{s,1}{params.channel_40x, 1});
    img_40x_mip = max(img_40x_mip, plane);
end
img_40x_mip = img_40x_mip ./ (max(img_40x_mip(:)) + eps);
[H40, W40]  = size(img_40x_mip);
fprintf('[40X]  MIP size: %d x %d px  (%.1f x %.1f um)  across %d Z planes\n', ...
        W40, H40, W40*px_40x, H40*px_40x, n_z);
%% -- Step 3: Resize 40X template to 10X pixel scale ----------------------
scale = px_40x / px_10x;   % 40X px / 10X px -- resize factor
fprintf('\n[REG]  Scale factor (40X->10X pixel space): %.4f\n', scale);
fprintf('[REG]  40X FOV physical size: %.1f x %.1f um\n', W40*px_40x, H40*px_40x);
% Resize MIP so it matches the 10X pixel scale
new_W = round(W40 * scale);
new_H = round(H40 * scale);
fprintf('[REG]  Resized template: %d x %d px in 10X space\n', new_W, new_H);
if new_W > W10 || new_H > H10
    error('[REG]  Resized 40X template (%dx%d) is larger than 10X image (%dx%d).\n       Check pixel sizes or series selection.', ...
          new_W, new_H, W10, H10);
end
template    = imresize(img_40x_mip, [new_H, new_W]);
template_eq = adapthisteq(template, 'NumTiles', [4 4], 'ClipLimit', 0.02);
%% -- Step 4: Test NCC against all 10X FOVs -------------------------------
fprintf('[REG]  Testing NCC against all %d 10X FOVs...\n', n_series_10x);
best_ncc    = -inf;
best_series = 1;
best_tl_row = 1;
best_tl_col = 1;
all_peaks   = zeros(1, n_series_10x);
all_xcorr   = cell(1, n_series_10x);
for s = 1:n_series_10x
    img_10x_s  = imgs_10x{s};
    img_10x_eq = adapthisteq(img_10x_s, 'NumTiles', [8 8], 'ClipLimit', 0.02);
    xcorr_s    = normxcorr2(template_eq, img_10x_eq);
    [peak_s, idx_s] = max(xcorr_s(:));
    all_peaks(s)    = peak_s;
    all_xcorr{s}    = xcorr_s;
    [pr, pc]        = ind2sub(size(xcorr_s), idx_s);
    fprintf('[REG]  FOV %2d: NCC peak = %.4f  (row=%d, col=%d)\n', s, peak_s, pr, pc);
    if peak_s > best_ncc
        best_ncc    = peak_s;
        best_series = s;
        best_tl_row = pr - new_H + 1;
        best_tl_col = pc - new_W + 1;
    end
end
fprintf('\n[REG]  Best match: 10X FOV %d  (NCC = %.4f)\n', best_series, best_ncc);
% Use best match for display
img_10x  = imgs_10x{best_series};
xcorr_map= all_xcorr{best_series};
tl_row   = best_tl_row;
tl_col   = best_tl_col;
br_row   = tl_row + new_H - 1;
br_col   = tl_col + new_W - 1;
[peak_val, peak_idx] = max(xcorr_map(:));
[peak_r, peak_c]     = ind2sub(size(xcorr_map), peak_idx);
cx_px = tl_col + new_W/2;
cy_px = tl_row + new_H/2;
cx_um = cx_px * px_10x;
cy_um = cy_px * px_10x;
fprintf('[REG]  FOV top-left  : (%d, %d) px\n',    tl_col, tl_row);
fprintf('[REG]  FOV centre    : (%.1f, %.1f) px\n', cx_px,  cy_px);
fprintf('[REG]  FOV centre    : (%.1f, %.1f) um\n', cx_um,  cy_um);
%% -- Step 5: Display results ----------------------------------------------
figure('Name', 'FOV Location Finder', 'NumberTitle', 'off', ...
       'Units', 'normalized', 'Position', [0.02 0.05 0.96 0.85]);
sgtitle(sprintf('40X FOV Location -- Best match: 10X FOV %d  (NCC=%.3f)', ...
        best_series, best_ncc), 'FontSize', 14, 'FontWeight', 'bold');
% Panel 1: Best 10X FOV with bounding box
ax1 = subplot(2,3,1);
imagesc(ax1, img_10x); colormap(ax1,'gray'); axis(ax1,'image'); axis(ax1,'off');
hold(ax1,'on');
r1c = max(1,tl_row); r2c = min(H10,br_row);
c1c = max(1,tl_col); c2c = min(W10,br_col);
rectangle('Position',[tl_col, tl_row, new_W, new_H], ...
          'EdgeColor',[0 1 0],'LineWidth',2);
plot(ax1, cx_px, cy_px, 'g+', 'MarkerSize',20, 'LineWidth',2);
title(ax1, sprintf('Best 10X FOV %d\nCentre: (%.0f, %.0f) um', ...
      best_series, cx_um, cy_um), 'FontSize',9);
% Panel 2: 40X MIP template
ax2 = subplot(2,3,2);
imagesc(ax2, img_40x_mip); colormap(ax2,'gray'); axis(ax2,'image'); axis(ax2,'off');
title(ax2, sprintf('40X MIP (%d Z planes)\n%.1f x %.1f um', ...
      n_z, W40*px_40x, H40*px_40x), 'FontSize',9);
% Panel 3: NCC scores across all 10X FOVs
ax3 = subplot(2,3,3);
bar(ax3, 1:n_series_10x, all_peaks, 'FaceColor', [0.3 0.6 1.0]);
hold(ax3,'on');
bar(ax3, best_series, all_peaks(best_series), 'FaceColor', [0.2 0.8 0.2]);
xlabel(ax3,'10X FOV index'); ylabel(ax3,'NCC peak score');
title(ax3, 'NCC score per 10X FOV (green=best)', 'FontSize',9);
grid(ax3,'on');
% Panel 4: NCC map for best match
ax4 = subplot(2,3,4);
imagesc(ax4, xcorr_map); colormap(ax4,'hot'); axis(ax4,'image'); axis(ax4,'off');
hold(ax4,'on');
plot(ax4, peak_c, peak_r, 'c+', 'MarkerSize',20,'LineWidth',2);
title(ax4, sprintf('NCC map -- FOV %d  (peak=%.3f)', best_series, peak_val),'FontSize',9);
colorbar(ax4);
% Panel 5: 10X crop at detected location
ax5 = subplot(2,3,5);
crop_10x = img_10x(r1c:r2c, c1c:c2c);
imagesc(ax5, crop_10x); colormap(ax5,'gray'); axis(ax5,'image'); axis(ax5,'off');
title(ax5, sprintf('10X crop -- FOV %d', best_series), 'FontSize',9);
% Panel 6: Side-by-side comparison
ax6 = subplot(2,3,6);
crop_resized = imresize(crop_10x, [new_H, new_W]);
comparison   = [mat2gray(crop_resized), mat2gray(template)];
imagesc(ax6, comparison); colormap(ax6,'gray'); axis(ax6,'image'); axis(ax6,'off');
title(ax6, '10X crop (left) vs 40X template (right)', 'FontSize',9);
drawnow;
%% -- Package output -------------------------------------------------------
result.best_series = best_series;
result.all_peaks   = all_peaks;
result.tl_col      = tl_col;
result.tl_row      = tl_row;
result.br_col      = br_col;
result.br_row      = br_row;
result.cx_px       = cx_px;
result.cy_px       = cy_px;
result.cx_um       = cx_um;
result.cy_um       = cy_um;
result.peak_ncc    = best_ncc;
result.scale       = scale;
result.px_10x      = px_10x;
result.px_40x      = px_40x;
result.xcorr_map   = xcorr_map;
fprintf('\n[DONE] Best match: 10X FOV %d  NCC=%.3f  centre=(%.1f, %.1f) um\n\n', ...
        best_series, best_ncc, cx_um, cy_um);
end
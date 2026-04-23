function display_results(results_L1, results_L2, results_L3, params)
%DISPLAY_RESULTS  Visualise results for Level 1, 2, 3, or any combination.
%
%  display_results(results_L1, results_L2, results_L3, params)
%  Pass [] for levels not run.
show_L1 = ~isempty(results_L1);
show_L2 = ~isempty(results_L2);
show_L3 = ~isempty(results_L3);
if show_L1; plot_single_level(results_L1, params); end
if show_L2; plot_single_level(results_L2, params); end
if show_L3; plot_level3(results_L3, params);       end
% Comparison: show whenever 2+ levels are run
active = {results_L1, results_L2, results_L3};
active = active(~cellfun(@isempty, active));
if numel(active) >= 2
    compare_levels(active, params);
end
end
%% ======================================================================
function plot_single_level(R, params)
%PLOT_SINGLE_LEVEL  Six-panel figure for Level 1 or Level 2.
fig = figure('Name', R.label, 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.90 0.85]);
sgtitle(sprintf('SAM Benchmark -- %s', R.label), ...
        'FontSize', 16, 'FontWeight', 'bold');
%% Panel 1 -- Segmentation mask
ax1 = subplot(2, 3, 1);
imagesc(ax1, R.seg_mask); colormap(ax1, 'jet');
axis(ax1, 'image'); axis(ax1, 'off');
title(ax1, sprintf('All Detected Nuclei\n(n=%d)', R.n_detected), 'FontSize', 10);
%% Panel 2 -- Scan pool
ax2 = subplot(2, 3, 2);
imagesc(ax2, R.seg_mask); colormap(ax2, 'gray');
axis(ax2, 'image'); axis(ax2, 'off'); hold(ax2, 'on');
for k = 1:numel(R.scan_cells)
    c = R.scan_cells(k).centroid;
    plot(ax2, c(1), c(2), '.', 'Color', [0 1 1], 'MarkerSize', 8);
end
title(ax2, sprintf('Scan Pool (n=%d)', R.n_scanned), 'FontSize', 10);
%% Panel 3 -- Selected targets
ax3 = subplot(2, 3, 3);
imagesc(ax3, R.seg_mask); colormap(ax3, 'gray');
axis(ax3, 'image'); axis(ax3, 'off'); hold(ax3, 'on');
for k = 1:numel(R.scan_cells)
    c = R.scan_cells(k).centroid;
    plot(ax3, c(1), c(2), '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5);
end
for k = 1:numel(R.target_cells)
    tc = R.target_cells(k);
    c  = tc.centroid;
    plot(ax3, c(1), c(2), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
    text(ax3, c(1)+8, c(2), sprintf('%.2f', tc.circularity), ...
         'Color', 'yellow', 'FontSize', 7, 'FontWeight', 'bold');
end
title(ax3, sprintf('Z-Stack Targets (n=%d)', R.n_targeted), 'FontSize', 10);
%% Panel 4 -- Circularity distribution
ax4 = subplot(2, 3, 4);
circ_all = [R.all_cells.circularity];
circ_sel = [R.target_cells.circularity];
h1 = histogram(ax4, circ_all, 20, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
hold(ax4, 'on');
h2 = histogram(ax4, circ_sel, 10, 'FaceColor', [1 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
xline(ax4, mean(circ_all), '--', 'Color', [0 0 0.8], 'LineWidth', 1.5, 'Label', 'mean all');
xline(ax4, mean(circ_sel), '-',  'Color', [0.8 0 0], 'LineWidth', 1.5, 'Label', 'mean sel.');
if R.level == 2
    xline(ax4, params.circularity_threshold, ':', 'Color', [0 0 0], 'LineWidth', 2, 'Label', 'threshold');
end
xlabel(ax4, 'Circularity'); ylabel(ax4, 'Count');
legend(ax4, [h1, h2], {'All cells', 'Selected'}, 'Location', 'best');
title(ax4, 'Circularity Distribution', 'FontSize', 10);
grid(ax4, 'on');
%% Panel 5 -- Z-stack strip
ax5 = subplot(2, 3, 5);
if ~isempty(R.zstack) && ~isempty(R.zstack(1).z_planes)
    zr = R.zstack(1); n_show = min(numel(zr.z_planes), 5);
    strip = [];
    for z = 1:n_show; strip = [strip, zr.z_planes{z}]; end %#ok<AGROW>
    imagesc(ax5, strip); colormap(ax5, 'gray');
    axis(ax5, 'image'); axis(ax5, 'off');
    title(ax5, sprintf('Cell %d -- Z-Stack Strip\n(best focus: Z%d)', ...
          zr.cell_id, zr.best_focus_z), 'FontSize', 10);
else
    axis(ax5, 'off');
    title(ax5, 'Z-Stack Preview', 'FontSize', 10);
end
%% Panel 6 -- Focus quality
ax6 = subplot(2, 3, 6);
if ~isempty(R.zstack) && ~isempty(R.zstack(1).focus_scores)
    hold(ax6, 'on');
    colors = lines(min(3, numel(R.zstack)));
    leg_h  = gobjects(min(3, numel(R.zstack)), 1);
    for k = 1:min(3, numel(R.zstack))
        zr = R.zstack(k);
        leg_h(k) = plot(ax6, 1:numel(zr.focus_scores), zr.focus_scores, ...
             '-o', 'Color', colors(k,:), 'LineWidth', 1.5, 'MarkerSize', 5);
        [~, bz] = max(zr.focus_scores);
        plot(ax6, bz, zr.focus_scores(bz), 'kv', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    end
    leg_labels = arrayfun(@(z) sprintf('Cell %d', z.cell_id), ...
                          R.zstack(1:min(3,end)), 'UniformOutput', false);
    legend(ax6, leg_h, leg_labels, 'Location', 'best', 'FontSize', 8);
    xlabel(ax6, 'Z Plane'); ylabel(ax6, 'Focus Score (Lap. var.)');
    title(ax6, 'Focus Quality per Z Plane', 'FontSize', 10);
    grid(ax6, 'on');
else
    axis(ax6, 'off');
end
drawnow;
save_figure(fig, sprintf('Level%d_%s.png', R.level, params.run_timestamp), params);
end
%% ======================================================================
function plot_level3(R, params)
%PLOT_LEVEL3  Six-panel figure for Level 3 (adds UNet confidence panel).
fig = figure('Name', R.label, 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.90 0.85]);
sgtitle(sprintf('SAM Benchmark -- %s', R.label), ...
        'FontSize', 16, 'FontWeight', 'bold');
%% Panel 1 -- Segmentation mask
ax1 = subplot(2, 3, 1);
imagesc(ax1, R.seg_mask); colormap(ax1, 'jet');
axis(ax1, 'image'); axis(ax1, 'off');
title(ax1, sprintf('Detected Nuclei (n=%d)', R.n_detected), 'FontSize', 10);
%% Panel 2 -- Scan pool coloured by UNet score
ax2 = subplot(2, 3, 2);
imagesc(ax2, R.seg_mask); colormap(ax2, 'gray');
axis(ax2, 'image'); axis(ax2, 'off'); hold(ax2, 'on');
scores = [R.scan_cells.unet_score];
cmap   = parula(256);
for k = 1:numel(R.scan_cells)
    c     = R.scan_cells(k).centroid;
    ci    = max(1, round(scores(k) * 255) + 1);
    color = cmap(ci, :);
    plot(ax2, c(1), c(2), 'o', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 8);
end
title(ax2, sprintf('Scan Pool -- coloured by UNet score\n(blue=low, yellow=high)'), 'FontSize', 10);
%% Panel 3 -- Selected targets
ax3 = subplot(2, 3, 3);
imagesc(ax3, R.seg_mask); colormap(ax3, 'gray');
axis(ax3, 'image'); axis(ax3, 'off'); hold(ax3, 'on');
for k = 1:numel(R.scan_cells)
    c = R.scan_cells(k).centroid;
    plot(ax3, c(1), c(2), '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 5);
end
for k = 1:numel(R.target_cells)
    tc = R.target_cells(k);
    c  = tc.centroid;
    plot(ax3, c(1), c(2), 'r+', 'MarkerSize', 14, 'LineWidth', 2);
    text(ax3, c(1)+8, c(2), sprintf('%.2f', tc.unet_score), ...
         'Color', 'yellow', 'FontSize', 7, 'FontWeight', 'bold');
end
xline_str = sprintf('threshold=%.2f', params.unet_threshold);
title(ax3, sprintf('Z-Stack Targets (n=%d)\nlabel=UNet score', R.n_targeted), 'FontSize', 10);
%% Panel 4 -- UNet confidence score distribution
ax4 = subplot(2, 3, 4);
all_scores = [R.scan_cells.unet_score];
sel_scores = [R.target_cells.unet_score];
h1 = histogram(ax4, all_scores, 15, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
hold(ax4, 'on');
h2 = histogram(ax4, sel_scores, 10, 'FaceColor', [1 0.4 0.1], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
xline(ax4, params.unet_threshold, ':', 'Color', [0 0 0], 'LineWidth', 2, 'Label', 'threshold');
xline(ax4, mean(sel_scores), '-', 'Color', [0.8 0 0], 'LineWidth', 1.5, 'Label', 'mean sel.');
xlabel(ax4, 'UNet Confidence Score'); ylabel(ax4, 'Count');
legend(ax4, [h1, h2], {'All scanned', 'Selected'}, 'Location', 'best');
title(ax4, sprintf('UNet Score Distribution\nModel: %s', params.unet_model), 'FontSize', 10);
grid(ax4, 'on');
%% Panel 5 -- Z-stack strip (high-res)
ax5 = subplot(2, 3, 5);
if ~isempty(R.zstack) && ~isempty(R.zstack(1).z_planes)
    zr = R.zstack(1); n_show = min(numel(zr.z_planes), 5);
    strip = [];
    for z = 1:n_show; strip = [strip, zr.z_planes{z}]; end %#ok<AGROW>
    imagesc(ax5, strip); colormap(ax5, 'gray');
    axis(ax5, 'image'); axis(ax5, 'off');
    title(ax5, sprintf('Cell %d -- High-Res Z-Stack\n(%d planes, best: Z%d)', ...
          zr.cell_id, numel(zr.z_planes), zr.best_focus_z), 'FontSize', 10);
else
    axis(ax5, 'off');
end
%% Panel 6 -- Focus quality
ax6 = subplot(2, 3, 6);
if ~isempty(R.zstack) && ~isempty(R.zstack(1).focus_scores)
    hold(ax6, 'on');
    colors = lines(min(3, numel(R.zstack)));
    leg_h  = gobjects(min(3, numel(R.zstack)), 1);
    for k = 1:min(3, numel(R.zstack))
        zr = R.zstack(k);
        leg_h(k) = plot(ax6, 1:numel(zr.focus_scores), zr.focus_scores, ...
             '-o', 'Color', colors(k,:), 'LineWidth', 1.5, 'MarkerSize', 5);
        [~, bz] = max(zr.focus_scores);
        plot(ax6, bz, zr.focus_scores(bz), 'kv', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    end
    leg_labels = arrayfun(@(z) sprintf('Cell %d', z.cell_id), ...
                          R.zstack(1:min(3,end)), 'UniformOutput', false);
    legend(ax6, leg_h, leg_labels, 'Location', 'best', 'FontSize', 8);
    xlabel(ax6, 'Z Plane'); ylabel(ax6, 'Focus Score');
    title(ax6, 'Focus Quality per Z Plane', 'FontSize', 10);
    grid(ax6, 'on');
else
    axis(ax6, 'off');
end
drawnow;
save_figure(fig, sprintf('Level3_%s.png', params.run_timestamp), params);
end
%% ======================================================================
function compare_levels(active_results, params)
%COMPARE_LEVELS  Side-by-side comparison for any combination of levels run.
n = numel(active_results);
fig = figure('Name', 'SAM Benchmark -- Level Comparison', 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.92 0.65]);
sgtitle('SAM Benchmark -- Level Comparison', 'FontSize', 15, 'FontWeight', 'bold');
colors = {[0.4 0.7 1.0], [1.0 0.4 0.4], [0.2 0.8 0.4]};
labels = cellfun(@(R) R.label, active_results, 'UniformOutput', false);
%% Row 1 -- Circularity distributions per level
for i = 1:n
    R  = active_results{i};
    ax = subplot(2, n, i);
    circ_all = [R.all_cells.circularity];
    circ_sel = [R.target_cells.circularity];
    h1 = histogram(ax, circ_all, 20, 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none');
    hold(ax, 'on');
    h2 = histogram(ax, circ_sel, 10, 'FaceColor', colors{i}, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    xline(ax, mean(circ_sel), '-', 'Color', colors{i}*0.6, 'LineWidth', 2);
    xlabel(ax, 'Circularity'); ylabel(ax, 'Count');
    title(ax, sprintf('%s\nmean sel. = %.3f', labels{i}, mean(circ_sel)), 'FontSize', 9);
    legend(ax, [h1 h2], {'All', 'Selected'}, 'Location', 'best', 'FontSize', 8);
    grid(ax, 'on');
end
%% Row 2 -- Bar chart of key metrics across all levels
ax_bar = subplot(2, n, (n+1):(2*n));
metric_labels = {'Circ. (selected)', 'Circ. (all)', 'Solidity', 'Aspect Ratio'};
vals = zeros(n, 4);
for i = 1:n
    R        = active_results{i};
    morph_i  = R.morph;
    vals(i,1) = mean([R.target_cells.circularity]);
    vals(i,2) = mean([R.all_cells.circularity]);
    vals(i,3) = mean([morph_i.solidity]);
    vals(i,4) = mean([morph_i.aspect_ratio]);
end
b = bar(ax_bar, vals', 'grouped');
for i = 1:n; b(i).FaceColor = colors{i}; end
set(ax_bar, 'XTickLabel', metric_labels, 'XTick', 1:4);
ylabel(ax_bar, 'Mean Value');
legend(ax_bar, labels, 'Location', 'best', 'FontSize', 9);
title(ax_bar, 'Morphology Metrics Comparison', 'FontSize', 11);
grid(ax_bar, 'on');
for i = 1:n
    for m = 1:4
        text(ax_bar, b(i).XEndPoints(m), b(i).YEndPoints(m)+0.01, ...
             sprintf('%.3f', vals(i,m)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
end
drawnow;
save_figure(fig, sprintf('Comparison_%s.png', params.run_timestamp), params);
%% Console summary
fprintf('\n+--------------------------------------------------------------+\n');
fprintf('| LEVEL COMPARISON SUMMARY                                    |\n');
fprintf('+----------------------+------------+------------+------------+\n');
fprintf('| Metric               |');
for i = 1:n; fprintf(' %-10s |', sprintf('L%d', active_results{i}.level)); end
fprintf('\n+----------------------+');
fprintf('%s\n', repmat('------------+', 1, n));
rows = {'Detected', 'Scanned', 'Targeted', 'Circ. sel.', 'Circ. all', 'Elapsed(s)'};
row_vals = zeros(n, numel(rows));
for i = 1:n
    R = active_results{i};
    row_vals(i,:) = [R.n_detected, R.n_scanned, R.n_targeted, ...
                     mean([R.target_cells.circularity]), ...
                     mean([R.all_cells.circularity]), R.elapsed_sec];
end
fmt = {'%d','%d','%d','%.3f','%.3f','%.2f'};
for r = 1:numel(rows)
    fprintf('| %-20s |', rows{r});
    for i = 1:n
        fprintf([' %-10' fmt{r}(end) ' |'], row_vals(i,r));
    end
    fprintf('\n');
end
fprintf('+----------------------+%s\n\n', repmat('------------+', 1, n));
end
%% ======================================================================
function save_figure(fig, filename, params)
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_path = fullfile(params.log_dir, filename);
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Figure saved: %s\n', filename);
end
end
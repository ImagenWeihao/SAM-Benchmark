function display_results(results_L1, results_L2, params)
%DISPLAY_RESULTS  Visualise benchmark results for Level 1, Level 2, or both.
%
%  display_results(results_L1, results_L2, params)
%  Pass [] for a level you did not run.

show_L1 = ~isempty(results_L1);
show_L2 = ~isempty(results_L2);

if show_L1; plot_single_level(results_L1, params); end
if show_L2; plot_single_level(results_L2, params); end
end


%% ======================================================================
function plot_single_level(R, params)
%PLOT_SINGLE_LEVEL  Six-panel figure for one benchmark level.

fig = figure('Name', R.label, 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.90 0.85]);
sgtitle(sprintf('SAM Benchmark -- %s', R.label), ...
        'FontSize', 16, 'FontWeight', 'bold');

%% Panel 1 -- Full FOV segmentation mask
ax1 = subplot(2, 3, 1);
imagesc(ax1, R.seg_mask);
colormap(ax1, 'jet');
axis(ax1, 'image'); axis(ax1, 'off');
title(ax1, sprintf('All Detected Nuclei\n(n=%d)', R.n_detected), 'FontSize', 10);

%% Panel 2 -- Scan pool
ax2 = subplot(2, 3, 2);
imagesc(ax2, R.seg_mask);
colormap(ax2, 'gray');
axis(ax2, 'image'); axis(ax2, 'off');
hold(ax2, 'on');
for k = 1:numel(R.scan_cells)
    c = R.scan_cells(k).centroid;
    plot(ax2, c(1), c(2), '.', 'Color', [0 1 1], 'MarkerSize', 8);
end
title(ax2, sprintf('Scan Pool\n(n=%d, cyan dots)', R.n_scanned), 'FontSize', 10);

%% Panel 3 -- Selected Z-stack targets
ax3 = subplot(2, 3, 3);
imagesc(ax3, R.seg_mask);
colormap(ax3, 'gray');
axis(ax3, 'image'); axis(ax3, 'off');
hold(ax3, 'on');
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
title(ax3, sprintf('Z-Stack Targets (n=%d)\nred=selected, label=circularity', ...
      R.n_targeted), 'FontSize', 10);

%% Panel 4 -- Circularity distribution
ax4 = subplot(2, 3, 4);
circ_all      = [R.all_cells.circularity];
circ_selected = [R.target_cells.circularity];

h1 = histogram(ax4, circ_all,      20, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
hold(ax4, 'on');
h2 = histogram(ax4, circ_selected, 10, 'FaceColor', [1 0.3 0.3],   'EdgeColor', 'none', 'FaceAlpha', 0.7);
xline(ax4, mean(circ_all),      '--', 'Color', [0 0 0.8], 'LineWidth', 1.5, 'Label', 'mean all');
xline(ax4, mean(circ_selected), '-',  'Color', [0.8 0 0], 'LineWidth', 1.5, 'Label', 'mean sel.');
if R.level == 2
    xline(ax4, params.circularity_threshold, ':', 'Color', [0 0 0], 'LineWidth', 2, 'Label', 'threshold');
end
xlabel(ax4, 'Circularity'); ylabel(ax4, 'Count');
legend(ax4, [h1, h2], {'All cells', 'Selected'}, 'Location', 'best');
title(ax4, 'Circularity Distribution', 'FontSize', 10);
grid(ax4, 'on');

%% Panel 5 -- Z-stack strip of first selected cell
ax5 = subplot(2, 3, 5);
if ~isempty(R.zstack) && ~isempty(R.zstack(1).z_planes)
    zr     = R.zstack(1);
    n_show = min(numel(zr.z_planes), 5);
    strip  = [];
    for z = 1:n_show
        strip = [strip, zr.z_planes{z}]; %#ok<AGROW>
    end
    imagesc(ax5, strip);
    colormap(ax5, 'gray');
    axis(ax5, 'image'); axis(ax5, 'off');
    title(ax5, sprintf('Cell %d -- Z-Stack Strip\n(planes 1 to %d, best focus: Z%d)', ...
          zr.cell_id, n_show, zr.best_focus_z), 'FontSize', 10);
else
    axis(ax5, 'off');
    text(0.5, 0.5, 'No Z-stack data', 'Units', 'normalized', ...
         'HorizontalAlignment', 'center', 'FontSize', 12);
    title(ax5, 'Z-Stack Preview', 'FontSize', 10);
end

%% Panel 6 -- Focus quality curves
ax6 = subplot(2, 3, 6);
if ~isempty(R.zstack) && ~isempty(R.zstack(1).focus_scores)
    hold(ax6, 'on');
    colors = lines(min(3, numel(R.zstack)));
    leg_handles = gobjects(min(3, numel(R.zstack)), 1);
    for k = 1:min(3, numel(R.zstack))
        zr = R.zstack(k);
        leg_handles(k) = plot(ax6, 1:numel(zr.focus_scores), zr.focus_scores, ...
             '-o', 'Color', colors(k,:), 'LineWidth', 1.5, 'MarkerSize', 5);
        [~, bz] = max(zr.focus_scores);
        plot(ax6, bz, zr.focus_scores(bz), 'kv', 'MarkerSize', 10, ...
             'MarkerFaceColor', 'k');
    end
    leg_labels = arrayfun(@(z) sprintf('Cell %d', z.cell_id), ...
                          R.zstack(1:min(3,end)), 'UniformOutput', false);
    legend(ax6, leg_handles, leg_labels, 'Location', 'best', 'FontSize', 8);
    xlabel(ax6, 'Z Plane'); ylabel(ax6, 'Focus Score (Lap. var.)');
    title(ax6, sprintf('Focus Quality per Z Plane\n(triangle = best focus)'), 'FontSize', 10);
    grid(ax6, 'on');
else
    axis(ax6, 'off');
    title(ax6, 'Focus Quality', 'FontSize', 10);
end

drawnow;

%% Save figure to log folder
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_name = sprintf('Level%d_%s.png', R.level, params.run_timestamp);
    fig_path = fullfile(params.log_dir, fig_name);
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Figure saved: %s\n', fig_path);
end
end


%% ======================================================================
function compare_levels(results_L1, results_L2, params)
%COMPARE_LEVELS  Side-by-side comparison figure: Level 1 vs Level 2.
% Called directly from run_benchmark when choice==3.

R1 = results_L1;
R2 = results_L2;

fig = figure('Name', 'SAM Benchmark -- L1 vs L2 Comparison', 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.05 0.90 0.70]);
sgtitle('Level 1 (Random) vs Level 2 (Circularity-Guided) Comparison', ...
        'FontSize', 15, 'FontWeight', 'bold');

levels = {R1, R2};
tstr   = {'L1: Random', sprintf('L2: Circ >= %.2f', params.circularity_threshold)};

%% Row 1 -- Circularity distributions per level
for i = 1:2
    R  = levels{i};
    ax = subplot(2, 2, i);
    circ_all = [R.all_cells.circularity];
    circ_sel = [R.target_cells.circularity];
    h1 = histogram(ax, circ_all, 20, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    hold(ax, 'on');
    h2 = histogram(ax, circ_sel, 10, 'FaceColor', [0.2 0.6 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    xline(ax, mean(circ_sel), '-', 'Color', [0.8 0 0], 'LineWidth', 2, 'Label', 'mean sel.');
    xlabel(ax, 'Circularity'); ylabel(ax, 'Count');
    title(ax, sprintf('%s\nSelected mean = %.3f', tstr{i}, mean(circ_sel)), 'FontSize', 10);
    legend(ax, [h1, h2], {'All', 'Selected'}, 'Location', 'best');
    grid(ax, 'on');
end

%% Row 2 -- Bar chart of key metrics
ax_bar = subplot(2, 2, 3:4);

morph1 = R1.morph;
morph2 = R2.morph;

metric_labels = {'Circ. selected', 'Circ. all', 'Solidity', 'Aspect Ratio'};
vals = [ mean([R1.target_cells.circularity]),  mean([R2.target_cells.circularity]);
         mean([R1.all_cells.circularity]),      mean([R2.all_cells.circularity]);
         mean([morph1.solidity]),               mean([morph2.solidity]);
         mean([morph1.aspect_ratio]),           mean([morph2.aspect_ratio]) ]';

b = bar(ax_bar, vals, 'grouped');
b(1).FaceColor = [0.4 0.7 1.0];
b(2).FaceColor = [1.0 0.4 0.4];
set(ax_bar, 'XTickLabel', metric_labels, 'XTick', 1:4);
ylabel(ax_bar, 'Mean Value');
legend(ax_bar, {'Level 1 (Random)', 'Level 2 (Circularity)'}, 'Location', 'best');
title(ax_bar, 'Morphology Metrics: L1 vs L2', 'FontSize', 11);
grid(ax_bar, 'on');

for g = 1:2
    for m = 1:4
        text(ax_bar, b(g).XEndPoints(m), b(g).YEndPoints(m) + 0.01, ...
             sprintf('%.3f', vals(g, m)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    end
end

drawnow;

%% Save comparison figure
if isfield(params, 'log_dir') && ~isempty(params.log_dir)
    fig_path = fullfile(params.log_dir, sprintf('Comparison_L1vsL2_%s.png', params.run_timestamp));
    exportgraphics(fig, fig_path, 'Resolution', 150);
    fprintf('[LOG]  Comparison figure saved: %s\n', fig_path);
end

%% Console summary table
fprintf('\n+----------------------------------------------------------------+\n');
fprintf('|          COMPARISON: Level 1 vs Level 2                       |\n');
fprintf('+--------------------------------+---------------+---------------+\n');
fprintf('| Metric                         | L1 (Random)   | L2 (Circ.)    |\n');
fprintf('+--------------------------------+---------------+---------------+\n');
fprintf('| Detected nuclei                | %-13d | %-13d |\n', R1.n_detected,  R2.n_detected);
fprintf('| Cells scanned                  | %-13d | %-13d |\n', R1.n_scanned,   R2.n_scanned);
fprintf('| Cells Z-stacked                | %-13d | %-13d |\n', R1.n_targeted,  R2.n_targeted);
fprintf('| Mean circ. (selected)          | %-13.3f | %-13.3f |\n', vals(1,1), vals(1,2));
fprintf('| Circ. improvement vs L1        | %-13s | +%-12.3f |\n', '--', vals(1,2)-vals(1,1));
fprintf('| Mean solidity                  | %-13.3f | %-13.3f |\n', vals(1,3), vals(2,3));
fprintf('| Mean aspect ratio              | %-13.3f | %-13.3f |\n', vals(1,4), vals(2,4));
fprintf('| Elapsed (s)                    | %-13.2f | %-13.2f |\n', R1.elapsed_sec, R2.elapsed_sec);
fprintf('+--------------------------------+---------------+---------------+\n\n');
end

function run_benchmark_gui()
%RUN_BENCHMARK_GUI  GUI front-end for the SAM benchmark.
%
%  Wraps run_level1 / run_level2 / run_level3 / run_segmentation_validation
%  with file+folder pickers and an optional real-time simulator mode that
%  drips .nd2 files in at a user-set interval.
%
%  Three configuration sections:
%    1. Input        : single .nd2 file OR folder of .nd2 files
%                      (+ simulator checkbox + interval spinner)
%    2. Ground Truth : enable checkbox + GT mask / CSV / nd2 paths
%                      + per-level L1/L2/L3 checkboxes
%    3. Results      : output folder (blank = auto logs\<timestamp>)
%
%  Launch:
%    >> run_benchmark_gui
%
%  The CLI workflow in run_benchmark.m is untouched.
%% -- Bio-Formats path setup (same as run_benchmark.m) --------------------
bf_path = 'D:\Matlab Files\toolbox\bfmatlab';
if isfolder(bf_path)
    addpath(bf_path);
    jar_file = fullfile(bf_path, 'bioformats_package.jar');
    existing_cp = javaclasspath('-dynamic');
    if ischar(existing_cp); existing_cp = {existing_cp}; end
    if exist(jar_file, 'file') && ~any(strcmp(existing_cp, jar_file))
        javaaddpath(jar_file);
    end
end
%% -- Initial state struct (captured in nested callbacks) -----------------
state.src_mode          = 'file';      % 'file' | 'folder'
state.src_path          = '';          % file path or folder path
state.simulator_on      = false;
state.sim_interval_s    = 5.0;
state.sim_timer         = [];
state.sim_queue         = {};          % remaining files in simulator run
state.gt_enabled        = false;
state.gt_mask           = '';
state.gt_csv            = '';
state.gt_nd2            = '';
state.run_l1            = true;
state.run_l2            = true;
state.run_l3            = true;
state.result_dir        = '';
state.running           = false;
%% -- Build the GUI -------------------------------------------------------
fig = figure('Name', 'SAM Benchmark Runner', ...
             'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none', ...
             'Position', [200 120 680 720], 'Resize', 'off', ...
             'CloseRequestFcn', @on_close);
h = struct();   % widget handles
% ---- Panel 1: Input --------------------------------------------------
p1 = uipanel(fig, 'Title', '1. Input .nd2', 'Units', 'pixels', ...
             'Position', [15 540 650 170], 'FontWeight', 'bold');
uicontrol(p1, 'Style', 'text', 'String', 'Source:', ...
          'Position', [10 130 60 20], 'HorizontalAlignment', 'left');
h.rb_file   = uicontrol(p1, 'Style', 'radiobutton', 'String', 'Single file', ...
                        'Value', 1, 'Position', [80 130 100 20], ...
                        'Callback', @(s,~) on_src_mode('file'));
h.rb_folder = uicontrol(p1, 'Style', 'radiobutton', 'String', 'Folder of .nd2', ...
                        'Value', 0, 'Position', [190 130 140 20], ...
                        'Callback', @(s,~) on_src_mode('folder'));
uicontrol(p1, 'Style', 'text', 'String', 'Path:', ...
          'Position', [10 95 40 20], 'HorizontalAlignment', 'left');
h.src_path = uicontrol(p1, 'Style', 'edit', 'String', '', ...
                       'Position', [50 95 490 24], ...
                       'HorizontalAlignment', 'left', 'Enable', 'off');
uicontrol(p1, 'Style', 'pushbutton', 'String', 'Browse...', ...
          'Position', [550 95 85 24], 'Callback', @on_browse_input);
h.cb_sim = uicontrol(p1, 'Style', 'checkbox', ...
          'String', 'Simulator mode  (drip files from folder at interval)', ...
          'Position', [10 60 400 20], 'Callback', @on_toggle_simulator);
uicontrol(p1, 'Style', 'text', 'String', 'Interval (s):', ...
          'Position', [10 30 70 20], 'HorizontalAlignment', 'left');
h.sim_interval = uicontrol(p1, 'Style', 'edit', 'String', '5.0', ...
                           'Position', [85 30 60 24], 'Enable', 'off');
% ---- Panel 2: Ground Truth ------------------------------------------
p2 = uipanel(fig, 'Title', '2. Ground Truth Benchmarking (optional)', ...
             'Units', 'pixels', 'Position', [15 310 650 220], ...
             'FontWeight', 'bold');
h.cb_gt = uicontrol(p2, 'Style', 'checkbox', ...
          'String', 'Enable GT benchmarking (validates against annotated mask)', ...
          'Position', [10 175 450 20], 'Callback', @on_toggle_gt);
uicontrol(p2, 'Style', 'text', 'String', 'GT mask (.tif):', ...
          'Position', [10 140 100 20], 'HorizontalAlignment', 'left');
h.gt_mask = uicontrol(p2, 'Style', 'edit', 'String', '', ...
                      'Position', [115 140 425 24], ...
                      'HorizontalAlignment', 'left', 'Enable', 'off');
h.btn_gt_mask = uicontrol(p2, 'Style', 'pushbutton', 'String', 'Browse...', ...
                          'Position', [550 140 85 24], 'Enable', 'off', ...
                          'Callback', @(s,~) on_browse_file('gt_mask', ...
                              {'*.tif;*.tiff', 'TIFF (*.tif)'}));
uicontrol(p2, 'Style', 'text', 'String', 'GT CSV:', ...
          'Position', [10 108 100 20], 'HorizontalAlignment', 'left');
h.gt_csv = uicontrol(p2, 'Style', 'edit', 'String', '', ...
                     'Position', [115 108 425 24], ...
                     'HorizontalAlignment', 'left', 'Enable', 'off');
h.btn_gt_csv = uicontrol(p2, 'Style', 'pushbutton', 'String', 'Browse...', ...
                         'Position', [550 108 85 24], 'Enable', 'off', ...
                         'Callback', @(s,~) on_browse_file('gt_csv', ...
                             {'*.csv', 'CSV (*.csv)'}));
uicontrol(p2, 'Style', 'text', 'String', 'GT .nd2:', ...
          'Position', [10 76 100 20], 'HorizontalAlignment', 'left');
h.gt_nd2 = uicontrol(p2, 'Style', 'edit', 'String', '', ...
                     'Position', [115 76 425 24], ...
                     'HorizontalAlignment', 'left', 'Enable', 'off');
h.btn_gt_nd2 = uicontrol(p2, 'Style', 'pushbutton', 'String', 'Browse...', ...
                         'Position', [550 76 85 24], 'Enable', 'off', ...
                         'Callback', @(s,~) on_browse_file('gt_nd2', ...
                             {'*.nd2', 'Nikon ND2 (*.nd2)'}));
uicontrol(p2, 'Style', 'text', 'String', 'Levels:', ...
          'Position', [10 40 60 20], 'HorizontalAlignment', 'left');
h.cb_l1 = uicontrol(p2, 'Style', 'checkbox', 'String', 'L1 (Otsu)', ...
                    'Position', [70 40 100 20], 'Value', 1, ...
                    'Callback', @(s,~) update_state('run_l1', logical(get(s,'Value'))));
h.cb_l2 = uicontrol(p2, 'Style', 'checkbox', 'String', 'L2 (Rule-based)', ...
                    'Position', [170 40 130 20], 'Value', 1, ...
                    'Callback', @(s,~) update_state('run_l2', logical(get(s,'Value'))));
h.cb_l3 = uicontrol(p2, 'Style', 'checkbox', 'String', 'L3 (StarDist)', ...
                    'Position', [300 40 120 20], 'Value', 1, ...
                    'Callback', @(s,~) update_state('run_l3', logical(get(s,'Value'))));
% ---- Panel 3: Results ------------------------------------------------
p3 = uipanel(fig, 'Title', '3. Result Folder', 'Units', 'pixels', ...
             'Position', [15 230 650 70], 'FontWeight', 'bold');
uicontrol(p3, 'Style', 'text', ...
          'String', 'Folder (blank = auto logs\<timestamp>):', ...
          'Position', [10 30 240 20], 'HorizontalAlignment', 'left');
h.result_dir = uicontrol(p3, 'Style', 'edit', 'String', '', ...
                         'Position', [255 28 285 24], ...
                         'HorizontalAlignment', 'left');
uicontrol(p3, 'Style', 'pushbutton', 'String', 'Browse...', ...
          'Position', [550 28 85 24], ...
          'Callback', @(s,~) on_browse_folder('result_dir'));
% ---- Control bar -----------------------------------------------------
h.btn_start = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Start', ...
                        'Position', [15 190 120 32], 'FontWeight', 'bold', ...
                        'BackgroundColor', [0.7 0.95 0.7], ...
                        'Callback', @on_start);
h.btn_stop  = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Stop', ...
                        'Position', [150 190 120 32], 'Enable', 'off', ...
                        'BackgroundColor', [0.95 0.75 0.7], ...
                        'Callback', @on_stop);
h.status    = uicontrol(fig, 'Style', 'text', 'String', 'Idle', ...
                        'Position', [290 195 370 24], ...
                        'HorizontalAlignment', 'left', ...
                        'BackgroundColor', [0.94 0.94 0.94]);
% ---- Log pane --------------------------------------------------------
uicontrol(fig, 'Style', 'text', 'String', 'Log:', ...
          'Position', [15 165 50 18], 'HorizontalAlignment', 'left', ...
          'FontWeight', 'bold');
h.log = uicontrol(fig, 'Style', 'listbox', 'String', {}, ...
                  'Position', [15 15 650 150], ...
                  'FontName', 'Consolas', 'FontSize', 9, ...
                  'BackgroundColor', [0.98 0.98 0.98]);
setappdata(fig, 'handles', h);
setappdata(fig, 'state', state);
log_msg(fig, 'GUI ready. Select input and click Start.');
%% ============================================================
%  Callbacks (all nested -- share state via appdata)
%% ============================================================
    function on_src_mode(mode)
        s = getappdata(fig, 'state');
        s.src_mode = mode;
        s.src_path = '';
        if strcmp(mode, 'file')
            set(h.rb_file,   'Value', 1); set(h.rb_folder, 'Value', 0);
            set(h.cb_sim, 'Enable', 'off', 'Value', 0);
            set(h.sim_interval, 'Enable', 'off');
            s.simulator_on = false;
        else
            set(h.rb_file,   'Value', 0); set(h.rb_folder, 'Value', 1);
            set(h.cb_sim, 'Enable', 'on');
        end
        set(h.src_path, 'String', '');
        setappdata(fig, 'state', s);
    end
    function on_browse_input(~, ~)
        s = getappdata(fig, 'state');
        if strcmp(s.src_mode, 'file')
            [file, folder] = uigetfile({'*.nd2', 'Nikon ND2 (*.nd2)'}, ...
                                       'Select input .nd2 file');
            if isequal(file, 0); return; end
            s.src_path = fullfile(folder, file);
        else
            folder = uigetdir('', 'Select folder containing .nd2 files');
            if isequal(folder, 0); return; end
            s.src_path = folder;
        end
        set(h.src_path, 'String', s.src_path);
        setappdata(fig, 'state', s);
        log_msg(fig, ['Input: ' s.src_path]);
    end
    function on_toggle_simulator(src, ~)
        s = getappdata(fig, 'state');
        s.simulator_on = logical(get(src, 'Value'));
        if s.simulator_on
            set(h.sim_interval, 'Enable', 'on');
        else
            set(h.sim_interval, 'Enable', 'off');
        end
        setappdata(fig, 'state', s);
    end
    function on_toggle_gt(src, ~)
        s = getappdata(fig, 'state');
        s.gt_enabled = logical(get(src, 'Value'));
        en = 'off'; if s.gt_enabled; en = 'on'; end
        set([h.gt_mask h.gt_csv h.gt_nd2 ...
             h.btn_gt_mask h.btn_gt_csv h.btn_gt_nd2], 'Enable', en);
        setappdata(fig, 'state', s);
    end
    function on_browse_file(field, filter)
        [file, folder] = uigetfile(filter, ['Select ' field]);
        if isequal(file, 0); return; end
        p = fullfile(folder, file);
        s = getappdata(fig, 'state');
        s.(field) = p;
        set(h.(field), 'String', p);
        setappdata(fig, 'state', s);
    end
    function on_browse_folder(field)
        folder = uigetdir('', ['Select ' field]);
        if isequal(folder, 0); return; end
        s = getappdata(fig, 'state');
        s.(field) = folder;
        set(h.(field), 'String', folder);
        setappdata(fig, 'state', s);
    end
    function update_state(field, val)
        s = getappdata(fig, 'state');
        s.(field) = val;
        setappdata(fig, 'state', s);
    end
    function on_start(~, ~)
        s = pull_gui_state(fig);
        if isempty(s.src_path)
            errordlg('Select an input file or folder first.', 'No input');
            return;
        end
        if s.gt_enabled && isempty(s.gt_mask)
            errordlg('GT benchmarking enabled but no GT mask selected.', ...
                     'Missing GT');
            return;
        end
        if ~s.run_l1 && ~s.run_l2 && ~s.run_l3
            errordlg('Enable at least one level (L1, L2, L3).', 'No level');
            return;
        end
        s.running = true;
        setappdata(fig, 'state', s);
        set(h.btn_start, 'Enable', 'off');
        set(h.btn_stop,  'Enable', 'on');
        if strcmp(s.src_mode, 'file')
            run_benchmark_on_file(fig, s.src_path);
            finish_run(fig);
        else
            start_folder_or_simulator(fig);
        end
    end
    function on_stop(~, ~)
        s = getappdata(fig, 'state');
        s.running = false;
        if ~isempty(s.sim_timer) && isvalid(s.sim_timer)
            stop(s.sim_timer);
            delete(s.sim_timer);
        end
        s.sim_timer = [];
        s.sim_queue = {};
        setappdata(fig, 'state', s);
        finish_run(fig);
        log_msg(fig, 'Stopped by user.');
    end
    function on_close(~, ~)
        s = getappdata(fig, 'state');
        if ~isempty(s.sim_timer) && isvalid(s.sim_timer)
            stop(s.sim_timer);
            delete(s.sim_timer);
        end
        delete(fig);
    end
end  % run_benchmark_gui
%% ============================================================
%  Helper functions (not nested, no closure access)
%% ============================================================
function s = pull_gui_state(fig)
%PULL_GUI_STATE  Refresh state from edit fields (user may have typed).
s = getappdata(fig, 'state');
h = getappdata(fig, 'handles');
s.src_path       = strtrim(get(h.src_path, 'String'));
s.sim_interval_s = max(0.1, str2double(get(h.sim_interval, 'String')));
if isnan(s.sim_interval_s); s.sim_interval_s = 5.0; end
s.gt_mask        = strtrim(get(h.gt_mask,    'String'));
s.gt_csv         = strtrim(get(h.gt_csv,     'String'));
s.gt_nd2         = strtrim(get(h.gt_nd2,     'String'));
s.result_dir     = strtrim(get(h.result_dir, 'String'));
setappdata(fig, 'state', s);
end
function start_folder_or_simulator(fig)
%Scan selected folder for .nd2 files, queue them, and run.
s = getappdata(fig, 'state');
files = dir(fullfile(s.src_path, '*.nd2'));
if isempty(files)
    log_msg(fig, sprintf('No .nd2 files found in %s', s.src_path));
    finish_run(fig);
    return;
end
queue = arrayfun(@(f) fullfile(f.folder, f.name), files, ...
                 'UniformOutput', false);
s.sim_queue = queue;
setappdata(fig, 'state', s);
log_msg(fig, sprintf('Queued %d .nd2 file(s).', numel(queue)));
if s.simulator_on
    log_msg(fig, sprintf('Simulator ON -- interval %.1fs', s.sim_interval_s));
    % Process first immediately, then drip rest on a timer.
    step_simulator(fig);
    t = timer('ExecutionMode', 'fixedSpacing', ...
              'Period', s.sim_interval_s, ...
              'BusyMode', 'queue', ...
              'TimerFcn', @(~,~) step_simulator(fig));
    s = getappdata(fig, 'state');
    s.sim_timer = t;
    setappdata(fig, 'state', s);
    start(t);
else
    % Process all sequentially (no delay between).
    while ~isempty(s.sim_queue)
        if ~s.running; break; end
        nd2_path = s.sim_queue{1};
        s.sim_queue(1) = [];
        setappdata(fig, 'state', s);
        run_benchmark_on_file(fig, nd2_path);
        drawnow;
        s = getappdata(fig, 'state');
    end
    finish_run(fig);
end
end
function step_simulator(fig)
%Called by timer: pops one file and runs it.
s = getappdata(fig, 'state');
if ~s.running || isempty(s.sim_queue)
    if ~isempty(s.sim_timer) && isvalid(s.sim_timer)
        stop(s.sim_timer);
        delete(s.sim_timer);
    end
    s.sim_timer = [];
    setappdata(fig, 'state', s);
    finish_run(fig);
    return;
end
nd2_path = s.sim_queue{1};
s.sim_queue(1) = [];
setappdata(fig, 'state', s);
log_msg(fig, sprintf('[SIM] New file arrives: %s', ...
        fileparts_name(nd2_path)));
run_benchmark_on_file(fig, nd2_path);
drawnow;
end
function run_benchmark_on_file(fig, nd2_path)
%Run the configured benchmark on a single .nd2 file.
s = getappdata(fig, 'state');
set_status(fig, sprintf('Processing %s ...', fileparts_name(nd2_path)));
log_msg(fig, sprintf('--- Start: %s ---', fileparts_name(nd2_path)));
try
    params = build_params(s, nd2_path);
    % Load once per file
    img_data = load_nd2(nd2_path, params);
    if ~s.gt_enabled
        % Plain level runs (L1/L2/L3)
        if s.run_l1
            log_msg(fig, 'Running L1 ...');
            run_level1(img_data, params);
        end
        if s.run_l2
            log_msg(fig, 'Running L2 ...');
            run_level2(img_data, params);
        end
        if s.run_l3
            log_msg(fig, 'Running L3 ...');
            run_level3(img_data, params);
        end
    else
        % GT validation route (equivalent to CLI Option 6)
        params.gt_mask_path = s.gt_mask;
        params.gt_csv_path  = s.gt_csv;
        params.gt_nd2_path  = s.gt_nd2;
        gt_params                   = params;
        gt_params.nucleus_channel   = params.gt_nucleus_channel;
        gt_params.series            = params.gt_series;
        gt_params.series_end        = params.gt_series_end;
        gt_params.z_plane           = params.gt_z_plane;
        gt_params.use_mip           = params.gt_use_mip;
        gt_img_data = load_nd2(s.gt_nd2, gt_params);
        log_msg(fig, 'Running GT validation ...');
        run_segmentation_validation(gt_img_data, s.gt_mask, s.gt_csv, params);
    end
    log_msg(fig, sprintf('Done: %s  ->  %s', ...
        fileparts_name(nd2_path), params.log_dir));
catch err
    log_msg(fig, sprintf('ERROR on %s: %s', ...
        fileparts_name(nd2_path), err.message));
end
set_status(fig, 'Idle');
end
function params = build_params(s, nd2_path)
%Construct the params struct the existing runners expect.
script_dir = fileparts(mfilename('fullpath'));
params = struct();
params.nd2_path        = nd2_path;
params.n_scan_cells    = 100;
params.n_zstack_cells  = 10;
params.nucleus_channel = 1;
params.gfp_channel     = 2;
params.rfp_channel     = 3;
params.far_red_channel = 4;
params.z_planes         = 7;
params.z_planes_highres = 15;
params.z_step_um        = 0.5;
params.circularity_threshold = 0.85;
% Level 3 / UNet
params.unet_model      = 'stardist';
params.unet_target     = 'round_nucleus';
params.unet_threshold  = 0.75;
params.unet_model_path = '';
% GT defaults (used when GT enabled)
params.gt_series          = 1;
params.gt_series_end      = 15;
params.gt_z_plane         = 1;
params.gt_nucleus_channel = 2;
params.gt_use_mip         = true;
params.gt_roi_enabled     = false;   % GUI path: skip interactive ROI
params.gt_export_centroids_csv = true;
% Result directory -- user-specified or auto
params.run_timestamp = datestr(now, 'yyyymmdd_HHMMSS');
[~, stem, ~] = fileparts(nd2_path);
if ~isempty(s.result_dir)
    params.log_dir = fullfile(s.result_dir, ...
        sprintf('%s_%s', stem, params.run_timestamp));
else
    params.log_dir = fullfile(script_dir, 'logs', ...
        sprintf('%s_%s', stem, params.run_timestamp));
end
if ~exist(params.log_dir, 'dir'); mkdir(params.log_dir); end
end
function finish_run(fig)
h = getappdata(fig, 'handles');
s = getappdata(fig, 'state');
s.running = false;
setappdata(fig, 'state', s);
set(h.btn_start, 'Enable', 'on');
set(h.btn_stop,  'Enable', 'off');
set(h.status, 'String', 'Idle');
end
function set_status(fig, msg)
h = getappdata(fig, 'handles');
set(h.status, 'String', msg);
drawnow;
end
function log_msg(fig, msg)
h = getappdata(fig, 'handles');
ts = datestr(now, 'HH:MM:SS');
entries = get(h.log, 'String');
if ischar(entries); entries = {entries}; end
entries{end+1} = sprintf('%s  %s', ts, msg);
if numel(entries) > 500; entries = entries(end-499:end); end
set(h.log, 'String', entries, 'Value', numel(entries));
drawnow;
end
function name = fileparts_name(p)
[~, n, ext] = fileparts(p);
name = [n ext];
end

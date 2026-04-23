function img_data = load_nd2(nd2_path, params)
%LOAD_ND2  Load a Nikon .nd2 file using Bio-Formats (bfopen).
%
%  img_data = load_nd2(nd2_path, params)
%
%  Returns a struct with fields:
%    .nucleus   - 2D grayscale image (nucleus channel)
%    .channels  - cell array {H x W double} one per channel
%    .metadata  - struct with pixel size, dimensions
%    .raw       - raw bfopen output
%% -- Check file exists ----------------------------------------------------
if ~isfile(nd2_path)
    error('SAM:fileNotFound', 'File not found:\n  %s', nd2_path);
end
%% -- Multi-file cache using containers.Map --------------------------------
% Stores one entry per unique nd2 file path -- loading two different files
% (e.g. main nd2 + GT nd2) does not evict each other from cache.
persistent cache_store
if isempty(cache_store)
    cache_store = containers.Map('KeyType','char','ValueType','any');
end
if isKey(cache_store, nd2_path)
    fprintf('[ND2]  Using cached data (skipping Bio-Formats reload).\n');
    data = cache_store(nd2_path);
else
    fprintf('[ND2]  Opening file with Bio-Formats...\n');
    data = bfopen(nd2_path);
    if isempty(data)
        error('SAM:bfEmpty', 'bfopen returned empty. File may be corrupt or unsupported.');
    end
    cache_store(nd2_path) = data;
    fprintf('[ND2]  File cached for future runs.\n');
end
%% -- Select series (FOV) --------------------------------------------------
if isfield(params, 'series') && params.series > 1
    s_idx = min(params.series, size(data, 1));
    fprintf('[ND2]  Using series %d of %d\n', s_idx, size(data, 1));
else
    s_idx = 1;
end
%% -- Extract metadata -----------------------------------------------------
omeMeta = data{1, 4};
try
    px_size_x = double(omeMeta.getPixelsPhysicalSizeX(s_idx-1).value());
    px_size_y = double(omeMeta.getPixelsPhysicalSizeY(s_idx-1).value());
catch
    try
        px_size_x = double(omeMeta.getPixelsPhysicalSizeX(0).value());
        px_size_y = double(omeMeta.getPixelsPhysicalSizeY(0).value());
    catch
        warning('Could not read pixel size. Assuming 0.325 um/px.');
        px_size_x = 0.325;
        px_size_y = 0.325;
    end
end
try
    n_channels = double(omeMeta.getPixelsSizeC(s_idx-1).getValue());
    n_z        = double(omeMeta.getPixelsSizeZ(s_idx-1).getValue());
    n_t        = double(omeMeta.getPixelsSizeT(s_idx-1).getValue());
    img_width  = double(omeMeta.getPixelsSizeX(s_idx-1).getValue());
    img_height = double(omeMeta.getPixelsSizeY(s_idx-1).getValue());
catch
    plane1     = data{s_idx,1}{1,1};
    img_height = size(plane1, 1);
    img_width  = size(plane1, 2);
    n_channels = params.far_red_channel;
    n_z        = max(1, floor(size(data{s_idx,1}, 1) / n_channels));
    n_t        = 1;
end
fprintf('[ND2]  Dimensions: W=%d  H=%d  C=%d  Z=%d  T=%d\n', ...
        img_width, img_height, n_channels, n_z, n_t);
fprintf('[ND2]  Pixel size: %.4f x %.4f um\n', px_size_x, px_size_y);
%% -- Extract per-channel images (single Z or MIP) -------------------------
image_stack = data{s_idx, 1};
channels    = cell(1, n_channels);
% MIP mode: compute Maximum Intensity Projection across multiple series
% Each series = one Z plane for this FOV (Nikon nd2 structure)
% Set params.use_mip = true and params.series_end to enable
use_mip    = isfield(params, 'use_mip') && params.use_mip;
series_end = s_idx;  % default: single series
if use_mip && isfield(params, 'series_end') && params.series_end > s_idx
    series_end = min(params.series_end, size(data, 1));
    fprintf('[ND2]  Computing MIP across series %d-%d (%d Z planes)...\n', ...
            s_idx, series_end, series_end - s_idx + 1);
elseif use_mip
    use_mip = false;   % only 1 series available, disable MIP
    fprintf('[ND2]  MIP requested but only 1 series -- loading single plane.\n');
end
if ~use_mip && isfield(params, 'z_plane') && params.z_plane > 1
    z_idx = min(params.z_plane, n_z);
    fprintf('[ND2]  Extracting Z plane %d of %d\n', z_idx, n_z);
else
    z_idx = 1;
end
for c = 1:n_channels
    if use_mip
        % Max-project channel c across all series in range
        mip_plane = zeros(img_height, img_width);
        for s = s_idx:series_end
            plane_idx = c;   % each series has 1 Z, so plane = channel index
            if plane_idx <= size(data{s,1}, 1)
                raw_z = double(data{s, 1}{plane_idx, 1});
                mip_plane = max(mip_plane, raw_z);
            end
        end
        raw_plane = mip_plane;
        z_idx     = 0;   % 0 = MIP
    else
        plane_idx = (z_idx - 1) * n_channels + c;
        if plane_idx <= size(image_stack, 1)
            raw_plane = image_stack{plane_idx, 1};
        else
            warning('Plane %d not found (Z=%d, C=%d); using zeros.', plane_idx, z_idx, c);
            raw_plane = zeros(img_height, img_width);
        end
    end
    raw_double = double(raw_plane);
    cmax = max(raw_double(:));
    if cmax > 0
        channels{c} = raw_double ./ cmax;
    else
        channels{c} = raw_double;
    end
end
%% -- Package output -------------------------------------------------------
img_data.nucleus    = channels{params.nucleus_channel};
img_data.channels   = channels;
img_data.n_channels = n_channels;
img_data.n_z        = n_z;
img_data.n_t        = n_t;
img_data.z_loaded   = z_idx;   % 0 = MIP, >0 = specific Z plane
img_data.raw        = data;
img_data.metadata   = struct( ...
    'px_size_x',  px_size_x, ...
    'px_size_y',  px_size_y, ...
    'img_width',  img_width, ...
    'img_height', img_height, ...
    'n_channels', n_channels, ...
    'n_z',        n_z, ...
    'n_t',        n_t, ...
    'z_loaded',   z_idx);
if use_mip
    fprintf('[ND2]  File loaded successfully (MIP across series %d-%d, %d Z planes).\n\n', ...
            s_idx, series_end, series_end - s_idx + 1);
else
    fprintf('[ND2]  File loaded successfully (Z=%d).\n\n', z_idx);
end
end
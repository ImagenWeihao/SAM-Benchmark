function img_data = load_nd2(nd2_path, params)
%LOAD_ND2  Load a Nikon .nd2 file using Bio-Formats (bfopen).
%
%  img_data = load_nd2(nd2_path, params)
%
%  Returns a struct with fields:
%    .nucleus   - 2D grayscale image (nucleus channel, first Z/T plane)
%    .channels  - cell array {H x W double} one per channel
%    .metadata  - struct with pixel size, dimensions, channel names
%    .raw       - raw bfopen output (for advanced use)

%% ── Check file exists ────────────────────────────────────────────────────
if ~isfile(nd2_path)
    error('SAM:fileNotFound', 'File not found:\n  %s', nd2_path);
end

%% ── Open file ────────────────────────────────────────────────────────────
fprintf('[ND2]  Opening file with Bio-Formats...\n');
data = bfopen(nd2_path);

if isempty(data)
    error('SAM:bfEmpty', 'bfopen returned empty. File may be corrupt or unsupported.');
end

%% ── Extract metadata ─────────────────────────────────────────────────────
omeMeta = data{1, 4};   % OME metadata object

try
    px_size_x = double(omeMeta.getPixelsPhysicalSizeX(0).value());
    px_size_y = double(omeMeta.getPixelsPhysicalSizeY(0).value());
catch
    warning('Could not read pixel size from metadata. Assuming 0.325 um/px.');
    px_size_x = 0.325;
    px_size_y = 0.325;
end

try
    n_channels = double(omeMeta.getPixelsSizeC(0).getValue());
    n_z        = double(omeMeta.getPixelsSizeZ(0).getValue());
    n_t        = double(omeMeta.getPixelsSizeT(0).getValue());
    img_width  = double(omeMeta.getPixelsSizeX(0).getValue());
    img_height = double(omeMeta.getPixelsSizeY(0).getValue());
catch
    plane1     = data{1,1}{1,1};
    img_height = size(plane1, 1);
    img_width  = size(plane1, 2);
    n_channels = params.far_red_channel;
    n_z        = max(1, floor(size(data{1,1}, 1) / n_channels));
    n_t        = 1;
end

fprintf('[ND2]  Dimensions: W=%d  H=%d  C=%d  Z=%d  T=%d\n', ...
        img_width, img_height, n_channels, n_z, n_t);
fprintf('[ND2]  Pixel size: %.4f x %.4f um\n', px_size_x, px_size_y);

%% ── Extract per-channel images (first Z plane, first T) ─────────────────
image_stack = data{1, 1};
channels    = cell(1, n_channels);

for c = 1:n_channels
    plane_idx = c;  % first Z, first T
    if plane_idx <= size(image_stack, 1)
        raw_plane = image_stack{plane_idx, 1};
    else
        warning('Plane %d not found; using zeros.', plane_idx);
        raw_plane = zeros(img_height, img_width);
    end
    raw_double = double(raw_plane);
    cmax = max(raw_double(:));
    if cmax > 0
        channels{c} = raw_double ./ cmax;
    else
        channels{c} = raw_double;
    end
end

%% ── Package output ───────────────────────────────────────────────────────
img_data.nucleus    = channels{params.nucleus_channel};
img_data.channels   = channels;
img_data.n_channels = n_channels;
img_data.n_z        = n_z;
img_data.n_t        = n_t;
img_data.raw        = data;
img_data.metadata   = struct( ...
    'px_size_x',  px_size_x, ...
    'px_size_y',  px_size_y, ...
    'img_width',  img_width, ...
    'img_height', img_height, ...
    'n_channels', n_channels, ...
    'n_z',        n_z, ...
    'n_t',        n_t);

fprintf('[ND2]  File loaded successfully.\n\n');
end

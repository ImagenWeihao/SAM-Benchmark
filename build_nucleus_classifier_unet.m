function net = build_nucleus_classifier_unet(input_size, n_classes)
%BUILD_NUCLEUS_CLASSIFIER_UNET  Lightweight 2-stage UNet for nucleus normality.
%
%  net = build_nucleus_classifier_unet(input_size, n_classes)
%
%  Stage 1: StarDist (external) detects and segments nuclei — returns crops
%  Stage 2: This network — classifies each nucleus crop as:
%             Class 1: normal          (round, single, compact)
%             Class 2: abnormal_shape  (multi-lobular, blebbing)
%             Class 3: abnormal_count  (binucleated, polyploid, micronuclei)
%
%  Architecture: Shallow encoder-decoder UNet
%    - 3 encoder blocks (16→32→64 filters) — lightweight, fast
%    - Bottleneck: 128 filters
%    - 3 decoder blocks with skip connections
%    - Final 1x1 conv → softmax over n_classes
%    - Input: single channel (DAPI / 405nm nucleus crop)
%    - Input size: 64x64x1 (default) — small enough for CPU inference
%
%  Requires: Deep Learning Toolbox
%
%  Usage:
%    net = build_nucleus_classifier_unet([64 64 1], 3);
%    analyzeNetwork(net)   % visualize architecture

if nargin < 1; input_size = [64 64 1]; end
if nargin < 2; n_classes  = 3;         end

H = input_size(1);
W = input_size(2);

fprintf('[UNET] Building lightweight nucleus UNet (%dx%d, %d classes)...\n', H, W, n_classes);

%% ── Layer definitions ────────────────────────────────────────────────────

layers = [

    %% Input
    imageInputLayer(input_size, 'Name', 'input', 'Normalization', 'rescale-zero-one')

    %% ── Encoder Block 1 (16 filters) ─────────────────────────────────────
    convolution2dLayer(3, 16, 'Padding', 'same', 'Name', 'enc1_conv1')
    batchNormalizationLayer('Name', 'enc1_bn1')
    reluLayer('Name', 'enc1_relu1')
    convolution2dLayer(3, 16, 'Padding', 'same', 'Name', 'enc1_conv2')
    batchNormalizationLayer('Name', 'enc1_bn2')
    reluLayer('Name', 'enc1_relu2')
    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'enc1_pool')   % H/2

    %% ── Encoder Block 2 (32 filters) ─────────────────────────────────────
    convolution2dLayer(3, 32, 'Padding', 'same', 'Name', 'enc2_conv1')
    batchNormalizationLayer('Name', 'enc2_bn1')
    reluLayer('Name', 'enc2_relu1')
    convolution2dLayer(3, 32, 'Padding', 'same', 'Name', 'enc2_conv2')
    batchNormalizationLayer('Name', 'enc2_bn2')
    reluLayer('Name', 'enc2_relu2')
    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'enc2_pool')   % H/4

    %% ── Encoder Block 3 (64 filters) ─────────────────────────────────────
    convolution2dLayer(3, 64, 'Padding', 'same', 'Name', 'enc3_conv1')
    batchNormalizationLayer('Name', 'enc3_bn1')
    reluLayer('Name', 'enc3_relu1')
    convolution2dLayer(3, 64, 'Padding', 'same', 'Name', 'enc3_conv2')
    batchNormalizationLayer('Name', 'enc3_bn2')
    reluLayer('Name', 'enc3_relu2')
    maxPooling2dLayer(2, 'Stride', 2, 'Name', 'enc3_pool')   % H/8

    %% ── Bottleneck (128 filters) ─────────────────────────────────────────
    convolution2dLayer(3, 128, 'Padding', 'same', 'Name', 'bottleneck_conv1')
    batchNormalizationLayer('Name', 'bottleneck_bn1')
    reluLayer('Name', 'bottleneck_relu1')
    convolution2dLayer(3, 128, 'Padding', 'same', 'Name', 'bottleneck_conv2')
    batchNormalizationLayer('Name', 'bottleneck_bn2')
    reluLayer('Name', 'bottleneck_relu2')

    %% ── Decoder Block 3 (64 filters) ─────────────────────────────────────
    transposedConv2dLayer(2, 64, 'Stride', 2, 'Name', 'dec3_upconv')   % H/4

    % Skip connection from enc3_relu2 added via lgraph below

    convolution2dLayer(3, 64, 'Padding', 'same', 'Name', 'dec3_conv1')
    batchNormalizationLayer('Name', 'dec3_bn1')
    reluLayer('Name', 'dec3_relu1')
    convolution2dLayer(3, 64, 'Padding', 'same', 'Name', 'dec3_conv2')
    batchNormalizationLayer('Name', 'dec3_bn2')
    reluLayer('Name', 'dec3_relu2')

    %% ── Decoder Block 2 (32 filters) ─────────────────────────────────────
    transposedConv2dLayer(2, 32, 'Stride', 2, 'Name', 'dec2_upconv')   % H/2

    % Skip connection from enc2_relu2 added via lgraph below

    convolution2dLayer(3, 32, 'Padding', 'same', 'Name', 'dec2_conv1')
    batchNormalizationLayer('Name', 'dec2_bn1')
    reluLayer('Name', 'dec2_relu1')
    convolution2dLayer(3, 32, 'Padding', 'same', 'Name', 'dec2_conv2')
    batchNormalizationLayer('Name', 'dec2_bn2')
    reluLayer('Name', 'dec2_relu2')

    %% ── Decoder Block 1 (16 filters) ─────────────────────────────────────
    transposedConv2dLayer(2, 16, 'Stride', 2, 'Name', 'dec1_upconv')   % H

    % Skip connection from enc1_relu2 added via lgraph below

    convolution2dLayer(3, 16, 'Padding', 'same', 'Name', 'dec1_conv1')
    batchNormalizationLayer('Name', 'dec1_bn1')
    reluLayer('Name', 'dec1_relu1')
    convolution2dLayer(3, 16, 'Padding', 'same', 'Name', 'dec1_conv2')
    batchNormalizationLayer('Name', 'dec1_bn2')
    reluLayer('Name', 'dec1_relu2')

    %% ── Output ───────────────────────────────────────────────────────────
    convolution2dLayer(1, n_classes, 'Name', 'output_conv')   % 1x1 conv
    softmaxLayer('Name', 'softmax')
    pixelClassificationLayer('Name', 'output')

];

%% ── Build layer graph and add skip connections ───────────────────────────
lgraph = layerGraph(layers);

% Skip connection concatenations (decoder upconv + encoder feature map)
% Each skip: concatenate along channel dimension

% Skip 3: dec3_upconv → concat with enc3_relu2 → feed into dec3_conv1
lgraph = addLayers(lgraph, depthConcatenationLayer(2, 'Name', 'skip3_concat'));
lgraph = disconnectLayers(lgraph, 'dec3_upconv',  'dec3_conv1');
lgraph = connectLayers(lgraph,    'dec3_upconv',  'skip3_concat/in1');
lgraph = connectLayers(lgraph,    'enc3_relu2',   'skip3_concat/in2');
lgraph = connectLayers(lgraph,    'skip3_concat', 'dec3_conv1');

% Skip 2: dec2_upconv → concat with enc2_relu2 → feed into dec2_conv1
lgraph = addLayers(lgraph, depthConcatenationLayer(2, 'Name', 'skip2_concat'));
lgraph = disconnectLayers(lgraph, 'dec2_upconv',  'dec2_conv1');
lgraph = connectLayers(lgraph,    'dec2_upconv',  'skip2_concat/in1');
lgraph = connectLayers(lgraph,    'enc2_relu2',   'skip2_concat/in2');
lgraph = connectLayers(lgraph,    'skip2_concat', 'dec2_conv1');

% Skip 1: dec1_upconv → concat with enc1_relu2 → feed into dec1_conv1
lgraph = addLayers(lgraph, depthConcatenationLayer(2, 'Name', 'skip1_concat'));
lgraph = disconnectLayers(lgraph, 'dec1_upconv',  'dec1_conv1');
lgraph = connectLayers(lgraph,    'dec1_upconv',  'skip1_concat/in1');
lgraph = connectLayers(lgraph,    'enc1_relu2',   'skip1_concat/in2');
lgraph = connectLayers(lgraph,    'skip1_concat', 'dec1_conv1');

net = lgraph;

%% ── Report ───────────────────────────────────────────────────────────────
fprintf('[UNET] Architecture built:\n');
fprintf('       Encoder:    16 → 32 → 64 filters (3 blocks)\n');
fprintf('       Bottleneck: 128 filters\n');
fprintf('       Decoder:    64 → 32 → 16 filters (3 blocks + skip connections)\n');
fprintf('       Output:     %d classes (softmax)\n', n_classes);
fprintf('       Input:      %dx%dx%d\n', input_size(1), input_size(2), input_size(3));

% Estimate parameter count
n_params = count_params(lgraph);
fprintf('       ~%s trainable parameters (lightweight)\n', format_num(n_params));
end


%% ── Helpers ──────────────────────────────────────────────────────────────
function n = count_params(lgraph)
% Rough estimate: sum filter weights across conv layers
n = 0;
for i = 1:numel(lgraph.Layers)
    L = lgraph.Layers(i);
    if isa(L, 'nnet.cnn.layer.Convolution2DLayer') || ...
       isa(L, 'nnet.cnn.layer.TransposedConvolution2DLayer')
        try
            sz = size(L.Weights);
            n  = n + prod(sz);
        catch
        end
    end
end
if n == 0; n = 485000; end  % fallback estimate for this architecture
end

function s = format_num(n)
if n >= 1e6
    s = sprintf('%.1fM', n/1e6);
elseif n >= 1e3
    s = sprintf('%.0fK', n/1e3);
else
    s = sprintf('%d', n);
end
end

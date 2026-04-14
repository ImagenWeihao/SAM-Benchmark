function [cells_out, model_info] = unet_mock(cells_in, img_data, params)
%UNET_MOCK  Deterministic mock classifier for testing without any ML model.
%
%  Scores cells using a combination of real morphology features:
%    circularity, intensity, and area — to produce a plausible confidence
%    score without requiring any model installation.
%
%  Useful for:
%    - Verifying the Level 3 pipeline runs end-to-end before installing models
%    - Baseline comparison against real ML models
%
%  Score formula:
%    score = 0.50 * circularity
%          + 0.30 * normalised_intensity
%          + 0.20 * area_score  (penalises size outliers)
%          + small Gaussian noise (seed=42, reproducible)

model_info.name    = 'Mock (morphology proxy)';
model_info.version = 'deterministic';
model_info.target  = params.unet_target;
model_info.status  = 'ok (mock)';

cells_out = cells_in;
n         = numel(cells_in);

% Normalise area to [0,1] within scan pool
areas      = [cells_in.area];
area_norm  = (areas - min(areas)) / (max(areas) - min(areas) + eps);
area_score = 1 - abs(area_norm - 0.5) * 2;   % penalise size outliers

% Normalise intensity
intens      = [cells_in.mean_int];
intens_norm = (intens - min(intens)) / (max(intens) - min(intens) + eps);

% Fixed seed for reproducibility
rng(42);
noise = 0.03 * randn(1, n);

for k = 1:n
    score = 0.50 * cells_in(k).circularity ...
          + 0.30 * intens_norm(k) ...
          + 0.20 * area_score(k) ...
          + noise(k);
    cells_out(k).unet_score = max(0, min(1, score));
    cells_out(k).unet_class = params.unet_target;
end

model_info.n_processed = n;
fprintf('[UNET] Mock scores computed from circularity + intensity + area.\n');
end

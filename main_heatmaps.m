clear, clc, close all
addpath('functions\');
addpath('sections\')
% snr_range      = 10 : 0.5 : 15;
% snr_range = [10, 13, 15];
snr_range = 10:0.5:15;
std_mult_range = 1  : 0.1 : 3.5;


nSNR  = length(snr_range);
nMult = length(std_mult_range);

% Pré-alocar matrizes de resultados [SNR × std_mult]
prec_map = zeros(nSNR, nMult);
rec_map  = zeros(nSNR, nMult);
f1_map   = zeros(nSNR, nMult);

conf_sweep = get_configuration(); 

[r, n] = rayleigh_profile(conf_sweep);

magnitudes = [+9.828e-6;   % +1.0 °C
              -4.914e-6;   % -0.5 °C
              +2.948e-6;   % +0.3 °C
              -9.828e-7;   % -0.1 °C
              +4.914e-6];  % +0.5 °C
% magnitudes = [-0.38e-6; 6.84e-7; -2.28e-7; +9.81e-7; -1.14e-7];
[gt, n_pert] = environmental_perturbation(conf_sweep, n, magnitudes);

fprintf('\nPulse shaping selection:     1) RC Filter Model      2) Super-Gaussian Model\n');
shape_choice = input('  > Option ');

if shape_choice == 1
    pulse_str = 'RC';
elseif shape_choice == 2
    pulse_str = 'SG';
else
    pulse_str = 'UnknownPulse';
end

timestamp_str = string(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));

fprintf('\nStarting SNR x std_mult sweep (%d x %d = %d simulations)...\n', ...
    nSNR, nMult, nSNR * nMult);

for i = 1:nSNR
    conf_sweep.SNR_dB = snr_range(i);

    % Recalcula backscatter só quando SNR muda
    [z_valid, t_laser, E_ref, E_sig, E_ref_raw_all] = ...
    backscatter_simulation(conf_sweep, shape_choice, r, n, n_pert);

    [freq_shift, freq_shift_calib, smooth_freq_shift, ~] = ...
    spectral_shift_estimation(conf_sweep, E_ref, E_sig, E_ref_raw_all);

    for j = 1:nMult
        conf_sweep.std_mult = std_mult_range(j);

        [all_locs, all_pks, ~] = peak_detection(freq_shift_calib, ...
            smooth_freq_shift, z_valid, conf_sweep);

        [Precision, Sensitivity, F1] = statistical_analysis(gt, conf_sweep, all_locs);

        prec_map(i,j) = Precision;
        rec_map (i,j) = Sensitivity;
        f1_map  (i,j) = F1;
    end

    fprintf('SNR = %.1f dB done!\n', snr_range(i));
end

%% --- Geração dos Heatmaps ---

% Criar uma nova figura

% 1. Heatmap da Precisão
fig8 = figure(8); set(gcf, 'Name', 'Precision Heatmap', 'WindowState', 'maximized');
h1 = heatmap(std_mult_range, snr_range, prec_map);
h1.Title = 'Precision';
h1.XLabel = 'std\_mult';
h1.YLabel = 'SNR (dB)';
h1.Colormap = jet; % Podes mudar para 'jet', 'hot', etc.
h1.ColorLimits = [0 100]; % Se estiver em %, ou [0 1] se forem frações
exportgraphics(fig8, fullfile('images', sprintf('%s_heatmap_precision_%s.png', ...
    pulse_str, timestamp_str)), 'Resolution', 300);

% 2. Heatmap do Recall (Sensitivity)
fig9 = figure(9); set(gcf, 'Name', 'Sensibility Heatmap', 'WindowState', 'maximized');
h2 = heatmap(std_mult_range, snr_range, rec_map);
h2.Title = 'Sensibility (Recall)';
h2.XLabel = 'std\_mult';
h2.YLabel = 'SNR (dB)';
h2.Colormap = jet;
h2.ColorLimits = [0 100];
exportgraphics(fig9, fullfile('images', sprintf('%s_heatmap_recall_%s.png', ...
    pulse_str, timestamp_str)), 'Resolution', 300);

% 3. Heatmap do F1-Score
fig10 = figure(10); set(gcf, 'Name', 'F1-Score Heatmap', 'WindowState', 'maximized');
h3 = heatmap(std_mult_range, snr_range, f1_map);
h3.Title = 'F1-Score';
h3.XLabel = 'std\_mult';
h3.YLabel = 'SNR (dB)';
h3.Colormap = jet; % Exemplo com outro mapa de cor para destacar o F1
h3.ColorLimits = [0 100];
exportgraphics(fig10, fullfile('images', sprintf('%s_heatmap_f1score_%s.png', ...
    pulse_str, timestamp_str)), 'Resolution', 300);
clear; clc; close all;

% Select and load the .mat data file
[file, path] = uigetfile('*.mat', 'Select the results file (.mat)');
if isequal(file, 0)
    disp('No file selected.');
    return;
end
load(fullfile(path, file));

%% FREQUENCY SHIFT
fig1 = figure(1);
set(fig1, 'Name', '1D Frequency Shift Profile');
z_axis = conf.z(1:conf.Nz - conf.M + 1);
plot(z_axis, freq_shift / 1e6, 'Color', [0.0 0.45 0.85], 'LineWidth', 1.2, 'DisplayName', 'Frequency Shift');
ax = gca;
set(ax, 'FontSize', 20, 'LineWidth', 0.8, 'FontName', 'Times New Roman');
xlabel('Distance (m)', 'FontSize', 20, 'FontName', 'Times New Roman'); 
ylabel('Frequency Shift (MHz)', 'FontSize', 20, 'FontName', 'Times New Roman');
%title('Distributed Frequency Shift Trace Along Sensing Fiber', ...
    %'FontSize', 22, 'FontWeight', 'bold');

xlim([207 237]);
ylim([-1100 1100]);
grid on;

%% PEAK DETECTIONS
% Separate positive and negative peaks
locs_pos = all_locs(all_pks > 0);
pks_pos  = all_pks(all_pks > 0);
locs_neg = all_locs(all_pks < 0);
pks_neg  = all_pks(all_pks < 0);

fig4 = figure(4);
set(fig4, 'Name', 'Peak Detection Diagnostics');
plot(z_valid, freq_shift/1e6, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'DisplayName', 'Raw trace');
hold on;
plot(z_valid, smooth_freq_shift/1e6, 'Color', [0.2 0.6 1.0], 'LineWidth', 2.5, 'DisplayName', 'Smoothed trace');
scatter(locs_pos, pks_pos/1e6, 250, 'g^', 'filled', 'DisplayName', 'Detected (+)');
scatter(locs_neg, pks_neg/1e6, 250, 'rv', 'filled', 'DisplayName', 'Detected (−)');
yline( threshold/1e6, 'g--', 'LineWidth', 1.8, 'DisplayName', '+Threshold');
yline(-threshold/1e6, 'r--', 'LineWidth', 1.8, 'DisplayName', '−Threshold');

% Positive Peaks
for p = 1:length(locs_pos)
    val_mhz = pks_pos(p)/1e6;
   
    y_pos = val_mhz + 80; % Offset to place label directly above peak
   
    txt = sprintf('%+.1f MHz', val_mhz);  % Label string
    
    text(locs_pos(p), y_pos, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'Color', [0.0 0.35 0.0], ...
        'BackgroundColor', [1 1 1 0.85], ...  % Semi-opaque white background
        'EdgeColor', [0.7 0.7 0.7], ...       % Subtle border box
        'Margin', 2, ...
        'FontName', 'Times New Roman');
end

% Negative Peaks
for n = 1:length(locs_neg)
    val_mhz = pks_neg(n)/1e6;
    
    y_pos = val_mhz - 80; % Offset to place label directly below peak
    
    txt = sprintf('%+.1f MHz', val_mhz);
    
    text(locs_neg(n), y_pos, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'Color', [0.5 0.0 0.0], ...
        'BackgroundColor', [1 1 1 0.85], ...
        'EdgeColor', [0.7 0.7 0.7], ...
        'Margin', 2,...
        'FontName', 'Times New Roman');
end
hold off;
set(gca, 'FontSize', 20, 'LineWidth', 1.2, 'Box', 'on', 'TickDir', 'in', 'FontName', 'Times New Roman');
xlabel('Distance (m)', 'FontSize', 20, 'FontName', 'Times New Roman'); 
ylabel('Frequency Shift (MHz)', 'FontSize', 20, 'FontName', 'Times New Roman');
% title('Distributed Spectral Shift Event Identification via Threshold Detection', ...
%     'FontSize', 22, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 14, 'FontName', 'Times New Roman'); 
grid on; 
xlim([207 240]); 
ylim([-900 1150]);

%% CROSS-CORRELATION
Nz = conf.Nz;
M = conf.M;
z = conf.z;
views = {
    [35, 45],   'Perspective View',    '35-45';
    [135, 45],  'Rear View',           '135-45';
    [0, 90],    'Top-Down / 2D View',  '0-90';
    [90, 40],   'Lateral View',        '90-40'
};

[Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);

for k = 1:size(views, 1)
    fig_view = figure(20 + k);
    set(fig_view, 'Name', sprintf('3D CC View %s', views{k,3}));
    surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none');
    colormap('jet');
    set(gca, 'FontSize', 25, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
    xlabel('Distance (m)', 'FontSize', 25, 'FontName', 'Times New Roman');
    ylabel('Frequency Lag (MHz)', 'FontSize', 25, 'FontName', 'Times New Roman');
    zlabel('Correlation', 'FontSize', 25);
    % title(sprintf('Cross-Correlation Surface Map — %s', views{k,2}), ...
        % 'FontSize', 30, 'FontWeight', 'bold');
    xlim([0 240]);
    ylim([-1100 1100]);
    zlim([-0.5 1]);
    view(views{k,1}(1), views{k,1}(2));
    rotate3d on;
    cb = colorbar;
    cb.FontSize = 20;
    
    % Dedicated adjustment for 90-40 view (prevents label overlap along Distance axis)
    set(gca, 'Position', [0.12, 0.15, 0.65, 0.72]); % Standard proportion
    cb.Position = [0.84, 0.15, 0.025, 0.72];
end

%%
Nz = conf.Nz;
M = conf.M;
z = conf.z;

figure(2)
[Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
contourf(Z_mesh, F_mesh, corr_map, 20, 'LineColor', 'none'); 
colormap('jet'); colorbar;
hold on;
% A linha branca traça o pico detetado pelo algoritmo matemático por cima do radar ótico
plot(z(1:Nz-M+1), freq_shift / 1e6, 'w', 'LineWidth', 1.5); 
hold off;
set(gca, 'FontSize', 25, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlabel('Distance (m)','FontSize', 25, 'FontName', 'Times New Roman'); 
ylabel('Frequency Lag (MHz)', 'FontSize', 25, 'FontName', 'Times New Roman');
% title('2D Correlation Map (Top View with Peak Trace)');
xlim([0 240]); 
ylim([-1000 1000]);

%%
Nz = conf.Nz;
M = conf.M;
z = conf.z;
z_valid = z(1 : (Nz - M + 1));

% =========================================================================
% FIGURE 5: IDEAL 3D MAP (NOISE-FREE)
% =========================================================================
fig5 = figure(5); set(fig5, 'Name', 'Cross-Correlation Ideal');
[Z_mesh, F_mesh] = meshgrid(z_valid(1:Nz-M+1), lags_freq / 1e6);
surf(Z_mesh, F_mesh, corr_map_id, 'EdgeColor', 'none');
view(35, 45); colormap('jet'); cb = colorbar; cb.FontSize = 12;
set(gca, 'FontSize', 15, 'LineWidth', 1.5);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Frequency Lag (MHz)', 'FontSize', 16);
zlabel('Correlation', 'FontSize', 16);
title('Ideal Spectral Cross-Correlation Map (Noise-Free Baseline)', ...
    'FontSize', 18, 'FontWeight', 'bold');
rotate3d on;
xlim([200 240]); ylim([-1100 1100]); zlim([-0.5 1]);
set(gca, 'Position', [0.12, 0.15, 0.65, 0.72]);
cb.Position = [0.84, 0.15, 0.025, 0.72];

%%
Nz = conf.Nz;
M = conf.M;
SNR_dB = conf.SNR_dB;
z = conf.z;
z_valid = z(1 : (Nz - M + 1));

% =========================================================================
% FIGURE 6: REAL 3D MAP (WITH NOISE)
% =========================================================================
fig6 = figure(6); set(fig6, 'Name', 'Cross-Correlation Real');
[Z_mesh, F_mesh] = meshgrid(z_valid(1:Nz-M+1), lags_freq / 1e6);
surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none');
view(35, 45); colormap('jet'); cb = colorbar; cb.FontSize = 12;
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'FontName', 'Times New Roman');
xlabel('Distance (m)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Frequency Lag (MHz)', 'FontSize', 20, 'FontName', 'Times New Roman');
zlabel('Correlation', 'FontSize', 20, 'FontName', 'Times New Roman');
% title(sprintf('Cross-Correlation Map under Noisy Conditions (SNR = %d dB)', ...
%     SNR_dB), 'FontSize', 18, 'FontWeight', 'bold');
rotate3d on;
xlim([200 240]); ylim([-1100 1100]); zlim([-0.5 1]);
set(gca, 'Position', [0.12, 0.15, 0.65, 0.72]); % Standard proportion
cb.Position = [0.84, 0.15, 0.025, 0.72];
cb.FontName = 'Times New Roman';
cb.FontSize = 20;

%%
Nz = conf.Nz;
M = conf.M;
z = conf.z;
z_valid = z(1 : (Nz - M + 1));

% =========================================================================
% FIGURE 7: RAW VS. SMOOTHED SIGNAL
% =========================================================================
smooth_window = M;
smooth_freq_shift = movmean(freq_shift, smooth_window);
fig7 = figure(7);
set(fig7, 'Name', 'Raw vs Smooth Frequency Shift');
plot(z_valid, freq_shift / 1e6, 'r-', ...
    'LineWidth', 1.5, 'DisplayName', 'Raw Trace', 'Color', [0.4 0.4 0.4]);
hold on;
plot(z_valid, smooth_freq_shift / 1e6, 'b-', ...
    'LineWidth', 1.5, 'DisplayName', 'Smoothed Trace', 'Color', [0.2 0.6 1.0]);
hold off;
grid on;
set(gca, 'FontSize', 15, 'LineWidth', 1.3);
xlabel('Distance (m)', 'FontSize', 16);
ylabel('Frequency Shift (MHz)', 'FontSize', 16);
title('Distributed Spectral Shift Profile (Raw vs. Smooth)', ...
    'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 14);
xlim([205, 240]);

%%
Nz = conf.Nz;
M = conf.M;
SNR_dB = conf.SNR_dB;
z = conf.z;
z_valid = z(1 : (Nz - M + 1));

% 1. Criar o vetor nulo para o Ground Truth contínuo
true_freq_shift_full = zeros(1, Nz);

% 2. Preencher com os desvios teóricos exatos de cada evento (em Hz)
for i = 1:conf.num_events
    start_idx = max(1, round(gt.theor_starts(i) / conf.dz));
    end_idx   = min(conf.Nz, round(gt.theor_ends(i) / conf.dz));
    true_freq_shift_full(start_idx:end_idx) = gt.theor_shifts(i) * 1e6; % Converte MHz para Hz
end

% 3. Truncar para alinhar com a dimensão de z_valid (Nz - M + 1)
true_freq_shift = true_freq_shift_full(1 : (conf.Nz - conf.M + 1));



% =========================================================================
% FIGURE 8: SMOOTHING WINDOW COMPARISON
% =========================================================================
smooth_window = M;
smooth_freq_shift = movmean(freq_shift, smooth_window);
smooth_small = movmean(freq_shift, max(1, round(M/2)));
smooth_large = movmean(freq_shift, 2*M);
fig8 = figure(8);
set(fig8, 'Name', 'Smoothing Window Trade-off Analysis');
z_gt = z_valid - (conf.M * conf.dz) / 2;

plot(z_valid, freq_shift / 1e6, '-', ...
    'LineWidth', 2.5, 'DisplayName', 'Raw Trace', 'Color', [0.7 0.7 0.7]);
hold on
plot(z_gt, true_freq_shift / 1e6, 'k--', ...
    'LineWidth', 2.5, 'DisplayName', 'Ground Truth');
plot(z_valid, smooth_small / 1e6, '--', ...
    'LineWidth', 2.5, 'DisplayName', sprintf('Under-smoothed (W = %d)', max(1, round(M/2))), ...
    'Color', [0.85 0.32 0.09]); % Orange
plot(z_valid, smooth_freq_shift / 1e6, '-.', ...
    'LineWidth',2.5, 'DisplayName', sprintf('Optimal (W = %d)', smooth_window), ...
    'Color', [0.00 0.45 0.74]); % Primary Blue
plot(z_valid, smooth_large / 1e6, '-', ...
    'LineWidth', 2.5, 'DisplayName', sprintf('Over-smoothed (W = %d)', 2*M), ...
    'Color', [0.47 0.67 0.19]); % Green
hold off;
grid on;
set(gca, 'FontSize', 20, 'LineWidth', 1.3, 'FontName', 'Times New Roman');
xlabel('Distance (m)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Frequency Shift (MHz)', 'FontSize', 20, 'FontName', 'Times New Roman');
%title('Trade-off of Moving-Average Window Size on Spectral Shift Profile', ...
    %'FontSize', 17, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 18, 'FontName', 'Times New Roman');
xlim([210, 226]);

%% HEATMAPS -- LOAD SECOND .MAT FILE
clear; clc; close all;

% Select and load the evaluation results .mat file
[file, path] = uigetfile('*.mat', 'Select the results file (.mat)');
if isequal(file, 0)
    disp('No file selected.');
    return;
end
load(fullfile(path, file));

%%
fig8 = figure(8); set(gcf, 'Name', 'Mean Precision Heatmap', 'WindowState', 'maximized');
h1 = heatmap(std_mult_range, snr_range, prec_map);
h1.Title = sprintf('Mean Precision (N = %d)', nExec);
h1.XLabel = 'Detection threshold scaling multiplier, k_{mult}';
h1.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h1.Colormap = jet;
h1.FontName = 'Times New Roman';
% h1.FontSize = 25;   22
h1.FontSize = 22;   
h1.ColorLimits = [0 100];
h1.CellLabelColor = 'none';
h1.GridVisible = 'off';
h1.YDisplayData = flip(h1.YDisplayData); % Inverte a ordem do Eixo Y

%%
fig9 = figure(9); set(gcf, 'Name', 'Mean Sensitivity Heatmap', 'WindowState', 'maximized');
h2 = heatmap(std_mult_range, snr_range, rec_map);
h2.Title = sprintf('Mean Sensitivity (N = %d)', nExec);
h2.XLabel = 'Detection threshold scaling multiplier, k_{mult}';
h2.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h2.Colormap = jet;
h2.FontSize = 22;
h2.ColorLimits = [0 100];
h2.FontName = 'Times New Roman';
h2.CellLabelColor = 'none';
h2.GridVisible = 'off';
h2.YDisplayData = flip(h2.YDisplayData);

%%
fig10 = figure(10); set(gcf, 'Name', 'Mean F1-Score Heatmap', 'WindowState', 'maximized');
h3 = heatmap(std_mult_range, snr_range, f1_map);
h3.Title = sprintf('Mean F1-Score (N = %d)', nExec);
h3.XLabel = 'Detection threshold scaling multiplier, k_{mult}';
h3.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h3.Colormap = jet;
h3.FontSize = 22;
h3.FontName = 'Times New Roman';
h3.ColorLimits = [0 100];
h3.CellLabelColor = 'none';
h3.GridVisible = 'off';
h3.YDisplayData = flip(h3.YDisplayData); % Inverte a ordem do Eixo Y

%%
fig11 = figure(11); set(gcf, 'Name', 'Precision Std Dev Heatmap', 'WindowState', 'maximized');
h4 = heatmap(std_mult_range, snr_range, prec_std_map);
h4.Title = sprintf('Precision Std. Dev. (N = %d)', nExec);
h4.XLabel = 'Detection threshold scaling multiplier, k_{mult}';
h4.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h4.Colormap = parula;
h4.FontSize = 22;
h4.FontName = 'Times New Roman';
h4.CellLabelColor = 'none';
h4.GridVisible = 'off';
h4.YDisplayData = flip(h4.YDisplayData); % Inverte a ordem do Eixo Y

%%
fig12 = figure(12); set(gcf, 'Name', 'Sensitivity Std Dev Heatmap', 'WindowState', 'maximized');
h5 = heatmap(std_mult_range, snr_range, rec_std_map);
h5.Title = sprintf('Sensitivity Std. Dev. (N = %d)', nExec);
h5.XLabel = 'Detection threshold scaling multiplier, k_{mult}';
h5.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h5.Colormap = parula;
h5.FontName = 'Times New Roman';
h5.FontSize = 22;
h5.CellLabelColor = 'none';
h5.GridVisible = 'off';
h5.YDisplayData = flip(h5.YDisplayData); % Inverte a ordem do Eixo Y

%%
fig13 = figure(13); set(gcf, 'Name', 'F1-Score Std Dev Heatmap', 'WindowState', 'maximized');
h6 = heatmap(std_mult_range, snr_range, f1_std_map);
h6.Title = sprintf('F1-Score Std. Dev. (N = %d)', nExec);
h6.XLabel = 'Detection threshold scaling multiplier, k_{mult}';
h6.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h6.FontName = 'Times New Roman';
h6.Colormap = parula;
h6.FontSize = 22;
h6.CellLabelColor = 'none';
h6.GridVisible = 'off';
h6.YDisplayData = flip(h6.YDisplayData); % Inverte a ordem do Eixo Y

%%
clear,close, clc

%%
addpath('functions\');
addpath('sections\');

% 1. Carregar configurações e perfil da fibra
conf = get_configuration();
[r, n] = rayleigh_profile(conf);

% Para esta figura, não precisamos de perturbação
n_pert = n; 

% 2. Executar a simulação de retroespalhamento
shape_choice = 3; % Impulso retangular ideal
[z_valid, ~, ~, ~, ~, ~, E_ref_raw_all, ~] = ...
    backscatter_simulation(conf, shape_choice, r, n, n_pert);

% 3. Extrair a intensidade |E(z)|^2 para uma única frequência ótica (ex: frequência central)
f_idx = round(conf.Nf / 2); % Índice da frequência central
E_z = E_ref_raw_all(f_idx, :); % Campo elétrico ao longo de z
intensity = abs(E_z).^2; % Intensidade |E(z)|^2

% NORMALIZAÇÃO DA INTENSIDADE (Atendendo à nota do orientador)
intensity_norm = intensity / max(intensity);

% 4. Gerar o gráfico do Traço de Rayleigh
fig50 = figure(50);
set(fig50, 'Name', 'Rayleigh Backscatter Trace');

% Plot da potência/intensidade normalizada
plot(z_valid, intensity_norm, 'LineWidth', 1.2);
grid on;

% Configuração dos eixos (Tamanho 20 + Fonte Times New Roman)
ax = gca;
set(ax, 'FontSize', 20, 'LineWidth', 1.3, 'FontName', 'Times New Roman');

xlabel('Distance (m)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Normalized Optical Intensity (a.u.)', 'FontSize', 20, 'FontName', 'Times New Roman');

xlim([0 conf.L]);
ylim([0 1.05]); % Define os limites verticais de 0 a 1 (com uma pequena margem no topo)

%%

%% HEATMAPS -- LOAD SECOND .MAT FILE
clear; clc; close all;

% Select and load the evaluation results .mat file
[file, path] = uigetfile('*.mat', 'Select the results file (.mat)');
if isequal(file, 0)
    disp('No file selected.');
    return;
end
load(fullfile(path, file));

%%
clc
% --- MÁXIMO ---
[max_val, max_idx] = max(f1_map(:));
[max_row, max_col] = ind2sub(size(f1_map), max_idx);

fprintf('==== VALOR MÁXIMO ====\n');
fprintf('F1-Score: %.2f%%\n', max_val);
fprintf('SNR: %g dB | k_mult: %g\n\n', snr_range(max_row), std_mult_range(max_col));

% --- MÍNIMO ---
[min_val, min_idx] = min(f1_map(:));
[min_row, min_col] = ind2sub(size(f1_map), min_idx);

fprintf('==== VALOR MÍNIMO ====\n');
fprintf('F1-Score: %.2f%%\n', min_val);
fprintf('SNR: %g dB | k_mult: %g\n', snr_range(min_row), std_mult_range(min_col));
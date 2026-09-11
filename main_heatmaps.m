%% Run N ITERATIONS — SNR x std_mult sweep with statistical robustness (Sec. 3.4.2)
clear, clc, close all
addpath('functions\');
addpath('sections\')

% ------------------------------------------------------------------
% GRID PARAMETERS
% ------------------------------------------------------------------
% Kept identical to the single-iteration sweep (Sec. 3.4.1) so that the
% mean-performance heatmaps here are directly comparable, cell-by-cell,
% to the single-run heatmaps. Reduce the grid resolution below (e.g.
% snr_range = 8:2:18, std_mult_range = 1.5:0.5:4.5) if runtime becomes
% impractical -- just be explicit about it in the thesis text if you do.
snr_range      = 8 : 1 : 18;      % dB, 11 points
std_mult_range = 1.5 : 0.2 : 4.5; % 16 points

nSNR  = length(snr_range);
nMult = length(std_mult_range);
nExec = 1;                       % iterations per (SNR, std_mult) pair
total_sims = nSNR * nMult * nExec;

% ------------------------------------------------------------------
% PRE-ALLOCATE RESULT MATRICES [SNR x std_mult]
% ------------------------------------------------------------------
prec_map = zeros(nSNR, nMult);
rec_map  = zeros(nSNR, nMult);
f1_map   = zeros(nSNR, nMult);

prec_std_map = zeros(nSNR, nMult);
rec_std_map  = zeros(nSNR, nMult);
f1_std_map   = zeros(nSNR, nMult);

conf_sweep = get_configuration();

[r, n] = rayleigh_profile(conf_sweep);

magnitudes = [+5.8968e-6;   % +0.6 
              -4.9140e-6;   % -0.5 
              +2.9484e-6;   % +0.3 
              -9.8280e-7;   % -0.1 
              +6.8796e-6];  % +0.7 

[gt, n_pert] = environmental_perturbation(conf_sweep, n, magnitudes);

fprintf(['\nPulse shaping selection:     \n' ...
    '   1) RC Filter Model\n' ...
    '   2) Super-Gaussian Model\n' ...
    '   3) Rectangular Pulse\n']);
shape_choice = input('> Option ');

if shape_choice == 1
    pulse_str = 'RC';
elseif shape_choice == 2
    pulse_str = 'SG';
elseif shape_choice == 3
    pulse_str = 'RP';
else
    pulse_str = 'UnknownPulse';
end

timestamp_str = string(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
out_dir = 'sg_window_31082026';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('\nStarting SNR x std_mult sweep (%d SNR x %d mult x %d exec = %d total simulations)...\n', ...
    nSNR, nMult, nExec, total_sims);

% ------------------------------------------------------------------
% MAIN SWEEP
% ------------------------------------------------------------------
main_tic = tic;
sim_count = 0;

for i = 1:nSNR
    conf_sweep.SNR_dB = snr_range(i);

    for j = 1:nMult
        conf_sweep.std_mult = std_mult_range(j);

        prec_runs = zeros(1, nExec);
        rec_runs  = zeros(1, nExec);
        f1_runs   = zeros(1, nExec);

        for k = 1:nExec
            sim_count = sim_count + 1;

            [z_valid, t_laser, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_all, E_ref_raw_calib_all] = ...
                backscatter_simulation(conf_sweep, shape_choice, r, n, n_pert);

            [freq_shift, freq_shift_calib, smooth_freq_shift, corr_map, lags_freq] = ...
                spectral_shift_estimation(conf_sweep, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_calib_all);

            [all_locs, all_pks, ~] = peak_detection(freq_shift_calib, ...
                smooth_freq_shift, z_valid, conf_sweep);

            [Precision, Sensitivity, F1] = statistical_analysis(gt, conf_sweep, all_locs);

            prec_runs(k) = Precision;
            rec_runs(k)  = Sensitivity;
            f1_runs(k)   = F1;
        end

        % Mean performance (comparable to the single-iteration maps)
        prec_map(i,j) = mean(prec_runs);
        rec_map(i,j)  = mean(rec_runs);
        f1_map(i,j)   = mean(f1_runs);

        % Run-to-run variability (new in the multi-iteration study)
        prec_std_map(i,j) = std(prec_runs);
        rec_std_map(i,j)  = std(rec_runs);
        f1_std_map(i,j)   = std(f1_runs);
    end

    elapsed_time = toc(main_tic);
    progress = sim_count / total_sims;
    avg_time_per_sim = elapsed_time / sim_count;
    remaining_time = avg_time_per_sim * (total_sims - sim_count);

    fprintf('SNR = %5.1f dB done! (Progress: %5.1f%% | Elapsed: %.1f min | ETA: %.1f min)\n', ...
        snr_range(i), progress * 100, elapsed_time/60, remaining_time/60);
end

total_duration = toc(main_tic);
fprintf('\n=========================================\n');
fprintf('Simulation Finished!\n');
fprintf('Total Execution Time: %.2f seconds (%.2f minutes)\n', total_duration, total_duration/60);
fprintf('Average time per simulation: %.3f seconds\n', total_duration/total_sims);
fprintf('=========================================\n');

%% ------------------------------------------------------------------
% SAVE RAW RESULTS (.mat) — keep before plotting in case of crash
% ------------------------------------------------------------------
mat_filename = fullfile(out_dir, sprintf('%s_sweep_results_N%d_%s.mat', ...
    pulse_str, nExec, timestamp_str));
save(mat_filename, 'prec_map', 'rec_map', 'f1_map', ...
    'prec_std_map', 'rec_std_map', 'f1_std_map', ...
    'snr_range', 'std_mult_range', 'nExec', 'total_duration');
fprintf('Raw results saved to: %s\n', mat_filename);

%% --- Mean performance heatmaps ---

fig8 = figure(8); set(gcf, 'Name', 'Mean Precision Heatmap', 'WindowState', 'maximized');
h1 = heatmap(std_mult_range, snr_range, prec_map);
h1.Title = sprintf('Mean Precision (N = %d)', nExec);
h1.XLabel = 'Detection threshold scaling multiplier, k\_mult';
h1.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h1.Colormap = jet;
h1.FontSize = 14;   
h1.ColorLimits = [0 100];
exportgraphics(fig8, fullfile(out_dir, sprintf('%s_heatmap_precision_mean_N%d_%s.png', ...
    pulse_str, nExec, timestamp_str)), 'Resolution', 300);

fig9 = figure(9); set(gcf, 'Name', 'Mean Sensitivity Heatmap', 'WindowState', 'maximized');
h2 = heatmap(std_mult_range, snr_range, rec_map);
h2.Title = sprintf('Mean Sensitivity / Recall (N = %d)', nExec);
h2.XLabel = 'Detection threshold scaling multiplier, k\_mult';
h2.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h2.Colormap = jet;
h2.FontSize = 14;
h2.ColorLimits = [0 100];
exportgraphics(fig9, fullfile(out_dir, sprintf('%s_heatmap_recall_mean_N%d_%s.png', ...
    pulse_str, nExec, timestamp_str)), 'Resolution', 300);

fig10 = figure(10); set(gcf, 'Name', 'Mean F1-Score Heatmap', 'WindowState', 'maximized');
h3 = heatmap(std_mult_range, snr_range, f1_map);
h3.Title = sprintf('Mean F1-Score (N = %d)', nExec);
h3.XLabel = 'Detection threshold scaling multiplier, k\_mult';
h3.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h3.Colormap = jet;
h3.FontSize = 14;
h3.ColorLimits = [0 100];
exportgraphics(fig10, fullfile(out_dir, sprintf('%s_heatmap_f1score_mean_N%d_%s.png', ...
    pulse_str, nExec, timestamp_str)), 'Resolution', 300);

%% --- Standard deviation (run-to-run variability) heatmaps ---

fig11 = figure(11); set(gcf, 'Name', 'Precision Std Dev Heatmap', 'WindowState', 'maximized');
h4 = heatmap(std_mult_range, snr_range, prec_std_map);
h4.Title = sprintf('Precision Std. Dev. (N = %d)', nExec);
h4.XLabel = 'Detection threshold scaling multiplier, k\_mult';
h4.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h4.Colormap = parula;
h4.FontSize = 14;
exportgraphics(fig11, fullfile(out_dir, sprintf('%s_heatmap_precision_std_N%d_%s.png', ...
    pulse_str, nExec, timestamp_str)), 'Resolution', 300);

fig12 = figure(12); set(gcf, 'Name', 'Sensitivity Std Dev Heatmap', 'WindowState', 'maximized');
h5 = heatmap(std_mult_range, snr_range, rec_std_map);
h5.Title = sprintf('Sensitivity Std. Dev. (N = %d)', nExec);
h5.XLabel = 'Detection threshold scaling multiplier, k\_mult';
h5.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h5.Colormap = parula;
h5.FontSize = 14;
exportgraphics(fig12, fullfile(out_dir, sprintf('%s_heatmap_recall_std_N%d_%s.png', ...
    pulse_str, nExec, timestamp_str)), 'Resolution', 300);

fig13 = figure(13); set(gcf, 'Name', 'F1-Score Std Dev Heatmap', 'WindowState', 'maximized');
h6 = heatmap(std_mult_range, snr_range, f1_std_map);
h6.Title = sprintf('F1-Score Std. Dev. (N = %d)', nExec);
h6.XLabel = 'Detection threshold scaling multiplier, k\_mult';
h6.YLabel = 'Signal-to-Noise Ratio, SNR (dB)';
h6.Colormap = parula;
h6.FontSize = 14;
exportgraphics(fig13, fullfile(out_dir, sprintf('%s_heatmap_f1score_std_N%d_%s.png', ...
    pulse_str, nExec, timestamp_str)), 'Resolution', 300);

fprintf('\nAll heatmaps exported to: %s\n', out_dir);
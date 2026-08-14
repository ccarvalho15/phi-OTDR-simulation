clear, clc, close all;
addpath('functions\');
addpath('sections\')

% 1. Inicia a gravação no ficheiro desejado
diary('report_rect_window_12082026.txt');

%% =======================================================================
% 1. SYSTEM CONFIGURATION & WAVEGU2IDE PROPERTIES
% ========================================================================
conf = get_configuration();

%% =======================================================================
% 2. STOCHASTIC RAYLEIGH SCATTERING PROFILE
% ========================================================================
[r, n] = rayleigh_profile(conf);

%% =======================================================================
% 3. ENVIRONMENTAL PERTURBATION MODEL (STRAIN / TEMPERATURE)
% ========================================================================

% Unique index change magnitudes per event 
%   positive = heating / tension,
%   negative = cooling / compression
% magnitudes = [-0.38e-6; 6.84e-7; -2.28e-7; +9.81e-7; -1.14e-7;];
magnitudes = [+5.8968e-6;   % +0.6 °C
              -4.9140e-6;   % -0.5 °C
              +2.9484e-6;   % +0.3 °C
              -9.8280e-7;   % -0.1 °C
              +6.8796e-6];  % +0.7 °C

[gt, n_pert] = environmental_perturbation(conf, n, magnitudes);

%% ========================================================================
% 4. & 5.                 PROBE SIGNAL & FREQUENCY SWEEP
%               COHERENT BACKSCATTER INTEGRATION & NOISE MODELLING
% ========================================================================

% PULSE SHAPE SELECTION
fprintf(['\nPulse shaping selection:     \n' ...
    '   1) RC Filter Model\n' ...
    '   2) Super-Gaussian Model\n' ...
    '   3) Rectangular Pulse\n']);
shape_choice = input('> Option ');

[z_valid, t_laser, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_all] = ...
    backscatter_simulation(conf, shape_choice, r, n, n_pert);
%% ========================================================================
% 6. SPECTRAL SHIFT ESTIMATION VIA CROSS-CORRELATION
% ========================================================================

[freq_shift, freq_shift_calib, smooth_freq_shift, corr_map, lags_freq] = ...
    spectral_shift_estimation(conf, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_all);


%% ========================================================================
% 7. PEAK DETECTION WITH NOISE SUPPRESSION
% ========================================================================

[all_locs, all_pks, threshold] = ...
    peak_detection(freq_shift_calib, smooth_freq_shift, z_valid, conf);

%% ========================================================================
% 8. PERFORMANCE REPORT
% ========================================================================

%   Calls the external 'report' utility which prints a confusion matrix,
%   per-event detection status, and figures-of-merit (precision, recall, F1).
report(gt.theor_starts, gt.theor_ends, gt.theor_shifts, gt.theor_delta_n, ...
    gt.theor_delta_T, all_locs, all_pks, z_valid, smooth_freq_shift, ...
    threshold, conf.gamma, conf.eta, conf.n_ave, conf.nu0, ...
    conf.pert_length, conf.match_tolerance)

%% ========================================================================
% 9. DATA VISUALISATION
% ========================================================================
%   This section visualizes the mapping between the physical perturbation and 
% the recovered frequency shifts, simulating the output of a distributed 
% fiber sensing interrogation system.

plot_results(conf, z_valid, freq_shift, smooth_freq_shift, ...
     corr_map, lags_freq, all_locs, all_pks, threshold);


%%
diary off;
fprintf('\n--- Simulation successfully completed! ---\n');


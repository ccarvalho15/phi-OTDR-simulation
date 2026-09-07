clear; clc; close all;
addpath('functions\');
addpath('sections\');

%% 1. CONFIGURATION & PARAMETERS
diary('rc_200_2026.09.06.txt');

conf = get_configuration();
num_iterations = 50;
num_events = conf.num_events;

% Pulse shape selection
fprintf(['\nPulse shaping selection:     \n' ...
    '   1) RC Filter Model\n' ...
    '   2) Super-Gaussian Model\n' ...
    '   3) Rectangular Pulse\n']);
shape_choice = input('> Option ');

% Spatial tolerance margin for peak assignment (meters)
match_tolerance = conf.match_tolerance;

% Refractive index perturbation magnitudes (ground truth)
magnitudes = [+5.8968e-6;   % +0.6 °C
              -4.9140e-6;   % -0.5 °C
              +2.9484e-6;   % +0.3 °C
              -9.8280e-7;   % -0.1 °C
              +6.8796e-6];  % +0.7 °C

% Spatial Localization Metrics Storage Matrices
mc_err_center = zeros(num_iterations, num_events);
mc_err_width  = zeros(num_iterations, num_events);
mc_err_w_rel  = zeros(num_iterations, num_events);

% Storage matrices per iteration [Iterations x Events]
mc_shifts_MHz = zeros(num_iterations, num_events);
mc_delta_n    = zeros(num_iterations, num_events);
mc_delta_T    = zeros(num_iterations, num_events);
mc_detected   = false(num_iterations, num_events);

fprintf('--- Starting Simulation (%d Iterations) ---\n\n', num_iterations);

%% 2. MONTE CARLO LOOP
main_tic = tic;
for iter = 1:num_iterations
    % a) Stochastic Rayleigh Profile and Environmental Model
    [r, n] = rayleigh_profile(conf);
    [gt, n_pert] = environmental_perturbation(conf, n, magnitudes);
    
    % b) Backscatter Signal Simulation
    [z_valid, t_laser, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_all, E_ref_raw_calib_all] = ...
        backscatter_simulation(conf, shape_choice, r, n, n_pert);
    
    % c) Spectral Shift Estimation (Cross-Correlation)
    [freq_shift, freq_shift_calib, smooth_freq_shift, corr_map, lags_freq] = ...
        spectral_shift_estimation(conf, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_calib_all);
    
    % d) Peak Detection via Threshold & Prominence
    [all_locs, all_pks, threshold] = peak_detection(freq_shift_calib, ...
        smooth_freq_shift, z_valid, conf);
    
    num_detected = length(all_locs);
    detected_matched = false(1, num_detected);
    
    % e) Spatial Event Matching & Classification Analysis
    for ev = 1:num_events
        theor_start  = gt.theor_starts(ev);
        theor_end    = gt.theor_ends(ev);
        theor_center = (theor_start + theor_end) / 2;
        
        limit_start = theor_start - match_tolerance;
        limit_end   = theor_end   + match_tolerance;
        
        % Search for detected peaks within spatial tolerance window
        match_idx = find(all_locs >= limit_start & all_locs <= limit_end);
        
        if ~isempty(match_idx)
            % Select the peak spatially closest to the event centroid
            [~, local_best] = min(abs(all_locs(match_idx) - theor_center));
            best_idx = match_idx(local_best);
            
            % True Positive (TP)
            detected_matched(best_idx) = true;
            mc_detected(iter, ev) = true;
            
            exp_peak_loc = all_locs(best_idx);
            
            % Trace spatial boundaries of the detected event
            [~, idx_peak] = min(abs(z_valid - exp_peak_loc));
            [idx_start, idx_end] = find_event_bounds(smooth_freq_shift, idx_peak, threshold);
            
            z_start = z_valid(idx_start);
            z_end   = z_valid(idx_end);
            exp_width = z_end - z_start;
            
            % Spatial Metrics Calculation
            mc_err_center(iter, ev) = exp_peak_loc - theor_center;
            mc_err_width(iter, ev)  = exp_width - conf.pert_length;
            mc_err_w_rel(iter, ev) = (mc_err_width(iter, ev) / conf.pert_length) * 100;
            
            % Extracted peak frequency shift (Hz to MHz)
            shift_Hz = all_pks(best_idx);
            mc_shifts_MHz(iter, ev) = shift_Hz / 1e6;
            
            % Physical Conversions: Delta_n and Delta_T
            dn_exp = (shift_Hz * conf.n_ave) / conf.nu0;
            dT_exp = temp_variation(dn_exp, conf.gamma, conf.eta, conf.n_ave);
            
            mc_delta_n(iter, ev) = dn_exp;
            mc_delta_T(iter, ev) = dT_exp;
        else
            mc_err_center(iter, ev) = NaN;
            mc_err_width(iter, ev)  = NaN;
            mc_err_w_rel(iter, ev)  = NaN;
            mc_shifts_MHz(iter, ev) = NaN;
            mc_delta_n(iter, ev)    = NaN;
            mc_delta_T(iter, ev)    = NaN;
        end
    end
    
    % --- Progress Tracking & ETA ---
    if rem(iter, 10) == 0 || iter == num_iterations
        elapsed = toc(main_tic);
        eta_min = (elapsed / iter) * (num_iterations - iter) / 60;
        fprintf('Iteration %3d/%3d completed (Progress: %5.1f%% | Elapsed: %.1f min | ETA: %.1f min)\n', ...
            iter, num_iterations, (iter/num_iterations)*100, elapsed/60, eta_min);
    end
end
total_duration = toc(main_tic);

%% 3. STATISTICAL METRICS CALCULATION
theor_starts     = gt.theor_starts;
theor_ends       = gt.theor_ends;
theor_shifts_MHz = gt.theor_shifts;
theor_delta_n    = gt.theor_delta_n;
theor_delta_T    = gt.theor_delta_T;

% Spatial Metrics Statistics
mean_err_center = mean(mc_err_center, 1, 'omitnan');
mean_err_width  = mean(mc_err_width, 1, 'omitnan');
mean_err_w_rel  = mean(mc_err_w_rel, 1, 'omitnan');

% Physical Parameters Statistics
mean_shifts = mean(mc_shifts_MHz, 1, 'omitnan');
std_shifts  = std(mc_shifts_MHz, 0, 1, 'omitnan');
err_shift_rel = abs(mean_shifts - theor_shifts_MHz) ./ abs(theor_shifts_MHz) * 100;

mean_dn = mean(mc_delta_n, 1, 'omitnan');
std_dn  = std(mc_delta_n, 0, 1, 'omitnan');
err_dn_rel = abs(mean_dn - theor_delta_n) ./ abs(theor_delta_n) * 100;

mean_dT = mean(mc_delta_T, 1, 'omitnan');
std_dT  = std(mc_delta_T, 0, 1, 'omitnan');
err_dT_rel = abs(mean_dT - theor_delta_T) ./ abs(theor_delta_T) * 100;

%% 4. DISPLAY COMPREHENSIVE RESULTS IN COMMAND WINDOW

% TABLE 1: LOCALIZATION (Theoretical vs Experimental positions)
fprintf('\n=========================================================================================\n');
fprintf('                                COMPARISON TABLE - LOCALIZATION\n');
fprintf('=========================================================================================\n');
fprintf('%-8s | %-16s | %-15s | %-14s | %-18s\n', ...
    'Event', 'Theoretical (m)', 'Err. Center (m)', 'Err. Width (m)', 'Rel. Error (%)');
fprintf('-----------------------------------------------------------------------------------------\n');
for ev = 1:num_events
    str_teo_loc = sprintf('[%.2f; %.2f]', theor_starts(ev), theor_ends(ev));
    fprintf('Event %-2d | %-16s | %-15.4f | %-14.4f | %-18.2f%%\n', ...
        ev, str_teo_loc, mean_err_center(ev), mean_err_width(ev), mean_err_w_rel(ev));
end

% TABLE 2: FREQUENCY SHIFT (MHz)
fprintf('\n=========================================================================================\n');
fprintf('                        COMPARISON TABLE - FREQUENCY SHIFT (MHz)\n');
fprintf('=========================================================================================\n');
fprintf('%-8s | %-17s | %-14s | %-13s | %-15s\n', 'Event', 'Theoretical (MHz)', 'Mean (MHz)', 'Std Dev', 'Rel. Error (%)');
fprintf('-----------------------------------------------------------------------------------------\n');
for ev = 1:num_events
    fprintf('Event %-2d | %+-17.2f | %+-14.2f | %-13.4f | %-15.2f%%\n', ...
        ev, theor_shifts_MHz(ev), mean_shifts(ev), std_shifts(ev), err_shift_rel(ev));
end

% TABLE 3: REFRACTIVE INDEX PERTURBATION (Delta n)
fprintf('\n=========================================================================================\n');
fprintf('                      COMPARISON TABLE - REFRACTIVE INDEX (DELTA N)\n');
fprintf('=========================================================================================\n');
fprintf('%-8s | %-17s | %-14s | %-13s | %-15s\n', 'Event', 'Theoretical', 'Mean', 'Std Dev', 'Rel. Error (%)');
fprintf('-----------------------------------------------------------------------------------------\n');
for ev = 1:num_events
    fprintf('Event %-2d | %+-17.4e | %+-14.4e | %-13.4e | %-15.2f%%\n', ...
        ev, theor_delta_n(ev), mean_dn(ev), std_dn(ev), err_dn_rel(ev));
end

% TABLE 4: TEMPERATURE VARIATION (Delta T in K)
fprintf('\n=========================================================================================\n');
fprintf('                      COMPARISON TABLE - TEMPERATURE VARIATION (DELTA T)\n');
fprintf('=========================================================================================\n');
fprintf('%-8s | %-17s | %-14s | %-13s | %-15s\n', 'Event', 'Theoretical (K)', 'Mean (K)', 'Std Dev', 'Rel. Error (%)');
fprintf('-----------------------------------------------------------------------------------------\n');
for ev = 1:num_events
    fprintf('Event %-2d | %+-17.4f | %+-14.4f | %-13.4f | %-15.2f%%\n', ...
        ev, theor_delta_T(ev), mean_dT(ev), std_dT(ev), err_dT_rel(ev));
end
fprintf('=========================================================================================\n');

%% 5. SAVE RESULTS TO MAT-FILE
output_dir = 'rc_window_2026.09.06_50s';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
filename  = fullfile(output_dir, sprintf('rc_200_report_%s.mat', timestamp));
save(filename, ...
    'mc_shifts_MHz', 'mc_delta_n', 'mc_delta_T', 'mc_detected', ...
    'mc_err_center', 'mc_err_width', 'mc_err_w_rel', ...
    'theor_shifts_MHz', 'theor_delta_n', 'theor_delta_T', ...
    'mean_err_center', 'mean_err_width', 'mean_err_w_rel', ...
    'mean_shifts', 'std_shifts', 'err_shift_rel', ...
    'mean_dn', 'std_dn', 'err_dn_rel', ...
    'mean_dT', 'std_dT', 'err_dT_rel', ...
    'gt', 'conf', 'num_iterations', 'total_duration');
fprintf('\nFull statistical report saved successfully to: %s\n', filename);
diary off;

%% ========================================================================
% HERLPER FUNCTIONS
% ========================================================================
function [idx_start, idx_end] = find_event_bounds(smooth_freq_shift, idx_peak, threshold)
% Dynamic edge-tracking algorithm via adaptive attenuation search.
%
%   Traces the boundaries of an event symmetrically outwards from its peak.
%   The boundaries are established at the spatial index where the energy of 
%   the smoothed frequency shift drops below 30% of the statistical 
%   threshold.
% -------------------------------------------------------------------------
    % Initialize boundaries at the localized peak index
    N = length(smooth_freq_shift);
    sign_peak = sign(smooth_freq_shift(idx_peak));
    idx_start = idx_peak;
    
    % Regress backward down the spatial trace until the frequency shift 
    % drops below 30% of the target noise floor threshold
    while idx_start > 1 && ...
            sign_peak * smooth_freq_shift(idx_start) > threshold
        idx_start = idx_start - 1;
    end
    
    % Progress forward down the spatial trace until the frequency shift
    % drops below 30% of the target noise floor threshold
    idx_end = idx_peak;
    while idx_end < N && ...
            sign_peak * smooth_freq_shift(idx_end) > threshold
        idx_end = idx_end + 1;
    end
end
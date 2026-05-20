function report(theor_starts, theor_ends, theor_shifts, theor_delta_n, ...
                               all_locs, all_pks, z_valid, smooth_freq_shift, ...
                               threshold, n_ave, nu0, pert_length, match_tolerance)
    % Get the total number of real/theoretical events and detected peaks
    num_events = length(theor_starts);
    num_detected = length(all_locs);
    
    %% ========================================================================
    % --- SECTION 7.5: ISOLATED DETECTION REPORT ---
    % ========================================================================
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('                DETECTION REPORT - PHI-OTDR\n');
    fprintf('%s\n', repmat('=', 1, 60));
    fprintf('%-8s | %-16s | %-17s | %-15s\n', 'Event', 'Location (m)', ...
        'Freq. Shift (MHz)', 'Delta_n');
    fprintf('%s\n', repmat('-', 1, 60));
    
    % Loop through each peak detected by the sensor to map its full physical width
    for i = 1:num_detected
        % Find the closest spatial index in the valid zone matching the detected peak location
        [~, idx_peak] = min(abs(z_valid - all_locs(i)));
        
        % Track backward and forward from the peak to find the start/end boundaries of the event.
        % The boundary is set where the frequency shift drops below 30% of the detection threshold.
        [idx_start, idx_end] = find_event_bounds(smooth_freq_shift, idx_peak, threshold);
        
        % Extract the physical fiber coordinates (in meters) for the calculated boundaries
        z_start = z_valid(idx_start);
        z_end   = z_valid(idx_end);
        loc_str = sprintf('[%.2f; %.2f]', z_start, z_end);
        
        % Invert the peak value to obtain the true physical frequency shift in Hz.
        % Sign convention: peak detection algorithms (like findpeaks) yield positive values,
        % but a physical compression/strain event produces a NEGATIVE optical frequency 
        % shift in phase-sensitive OTDR. Negation recovers the true signed shift.
        shift_Hz = - all_pks(i); 
        
        % Calculate the refractive index variation (Delta_n) using the optical relationship:
        % Delta_n = (Delta_nu * n_average) / nu_0
        delta_n_val = (shift_Hz * n_ave) / nu0;
        
        % Print individual detected event parameters (converting Hz to MHz for readability)
        fprintf('Event %-2d | %-16s | %+-17.2f | %+-15.2e\n', ...
            i, loc_str, shift_Hz/1e6, delta_n_val);
    end
    fprintf('%s\n', repmat('-', 1, 60));
    
    %% ========================================================================
    % --- STATISTICAL PERFORMANCE ANALYSIS (TP, FP, FN) ---
    % ========================================================================
    detected_matched = false(1, num_detected); % Track which detected peaks match actual ground truth
    TP = 0;                        % True Positives counter (correctly identified events)
    FN = 0;                        % False Negatives counter (missed events)
    
    fprintf('\n\n%s\n', repmat('=', 1, 66));
    fprintf('         STATISTICAL PERFORMANCE ANALYSIS (TP, FP, FN)\n');
    fprintf('%s\n', repmat('=', 1, 66));
    
    % Evaluate True Positives (TP) and False Negatives (FN) against ground truth
    for i = 1:num_events
        % Define spatial acceptance boundaries using the tolerance margin
        limit_start = theor_starts(i) - match_tolerance;
        limit_end = theor_ends(i) + match_tolerance;
        
        % Check if any detected peaks fall within the valid spatial window
        match_idx = find(all_locs >= limit_start & all_locs <= limit_end);
        
        if ~isempty(match_idx)
            TP = TP + 1;
            detected_matched(match_idx) = true; % Flag these peaks as successfully matched
            fprintf('TP (Event %d at [%.2f; %.2f] m) :: DETECTED at peak %.2f m\n', ...
                i, theor_starts(i), theor_ends(i), all_locs(match_idx(1)));
        else
            FN = FN + 1;
            fprintf('FN (Event %d at [%.2f; %.2f] m) :: NOT DETECTED by the sensor\n', ...
                i, theor_starts(i), theor_ends(i));
        end
    end
    
    % Evaluate False Positives (FP)
    % False Positives are detected peaks that did not map to any real ground truth event
    FP = sum(~detected_matched);
    if FP > 0
        fp_indices = find(~detected_matched);
        for k = 1:length(fp_indices)
            fprintf('FP (Noise Artifact) :: DETECTED at peak %.2f m\n', all_locs(fp_indices(k)));
        end
    else
        fprintf('No False Alarms (FP) detected outside the perturbation zones.\n');
    end
    
    % Compute percentages for the confusion matrix parameters
    pct_TP = (TP / num_events) * 100;
    pct_FN = (FN / num_events) * 100;
    if num_detected > 0
        pct_FP = (FP / num_detected) * 100;
    else
        pct_FP = 0;
    end
    
    % Calculate key signal detection metrics: Sensitivity, Precision, and F1-Score
    Sensitivity = (TP / (TP + FN)) * 100; 
    if (TP + FP) > 0
        Precision = (TP / (TP + FP)) * 100;
    end
    if (Precision + Sensitivity) > 0
        F1 = 2 * (Precision * Sensitivity) / (Precision + Sensitivity);
    end
    
    % Print the comprehensive classification performance report
    fprintf('%s\n', repmat('=', 1, 66));
    fprintf('                       SUMMARIZE REPORT\n');
    fprintf('%s\n', repmat('=', 1, 66));
    fprintf('Total Real Perturbations: %d\n', num_events);
    fprintf('Total Detected Peaks:     %d\n', num_detected);
    fprintf('%s\n', repmat('-', 1, 66));
    fprintf('True Positives (TP):      %d (%.2f%% of real events detected)\n', TP, pct_TP);
    fprintf('False Negatives (FN):     %d (%.2f%% of real events missed)\n', FN, pct_FN);
    fprintf('False Positives (FP):     %d (%.2f%% of triggers were false alarms)\n', FP, pct_FP);
    fprintf('%s\n', repmat('-', 1, 66));
    fprintf('System Sensitivity:       %.1f%%\n', Sensitivity);
    fprintf('System Precision:         %.1f%%\n', Precision);
    fprintf('F1-Score:                 %.1f%%\n', F1);
    fprintf('%s\n\n\n', repmat('=', 1, 66));
    
    %% ========================================================================
    % --- SECTION 7.6: DIRECT COMPARATIVE TRANSLATION PROCESSING ---
    % ========================================================================
    % Pre-allocate memory arrays for error tracking (Experimental vs. Theoretical)
        
    % Spatial Localization Arrays
    err_center = NaN(num_events, 1);     % Centroid location error
    err_width = NaN(num_events, 1);      % Broadening/narrowing error (width mismatch)
    err_width_rel = NaN(num_events, 1);  % Relative width error (%)
    str_det_loc = cell(num_events, 1);   % Formatted text storage for [start; end] boundaries
    
    % Optical Frequency Shift Arrays
    shift_MHz_exp = zeros(num_events, 1); % Captured experimental shift in MHz
    err_shift_abs = NaN(num_events, 1);   % Absolute error of frequency shift
    err_shift_rel = NaN(num_events, 1);   % Relative error of frequency shift (%)
    
    % Refractive Index Mismatch Arrays
    dn_exp = zeros(num_events, 1);       % Calculated experimental Delta_n
    err_dn_abs = NaN(num_events, 1);     % Absolute error of Delta_n
    err_dn_rel = NaN(num_events, 1);     % Relative error of Delta_n (%)
    
    % Detection status flags
    is_detected = false(num_events, 1);
    
    % Evaluate errors by performing one-to-one comparisons for each real event
    for i = 1:num_events
        % Isolate the index range of the valid fibers matching the current theoretical zone
        idx_theor_zone = find(z_valid >= theor_starts(i) & z_valid <= theor_ends(i));
        
        theor_peak_lock = NaN; 
        if ~isempty(idx_theor_zone)
            % Locate the point of maximum absolute frequency shift inside this specific ground-truth window
            [~, local_max_idx] = max(abs(smooth_freq_shift(idx_theor_zone)));
            theor_peak_lock = z_valid(idx_theor_zone(local_max_idx)); % Pinpoint theoretical peak location
        end
        
        if ~isempty(all_locs) && ~isnan(theor_peak_lock)
            % Identify which detected peak is spatially closest to the targeted theoretical peak
            [dist_min, best_idx] = min(abs(all_locs - theor_peak_lock));
            
            % Association is verified only if the closest peak falls within the allowed spatial tolerance radius
            if dist_min < match_tolerance
                is_detected(i) = true;
                exp_peak_loc = all_locs(best_idx);
                
                % Dynamically evaluate event edge boundaries to gauge broadening errors
                [~, idx_peak] = min(abs(z_valid - exp_peak_loc));
                [idx_start, idx_end] = find_event_bounds(smooth_freq_shift, idx_peak, threshold);
                
                z_start = z_valid(idx_start);
                z_end   = z_valid(idx_end);
                str_det_loc{i} = sprintf('[%.2f; %.2f]', z_start, z_end);
                
                % Calculate spatial error metrics (difference between peak centers and event widths)
                err_center(i) = exp_peak_loc - theor_peak_lock;
                exp_width     = z_end - z_start;        
                err_width(i)  = exp_width - pert_length;
                if pert_length ~= 0
                    err_width_rel(i) = (err_width(i) / pert_length) * 100;
                end
                
                % Calculate frequency shift error metrics (converting experimental data to MHz)
                shift_MHz_exp(i) = -all_pks(best_idx) / 1e6; 
                err_shift_abs(i) = abs(shift_MHz_exp(i) - theor_shifts(i));
                if theor_shifts(i) ~= 0
                    err_shift_rel(i) = (err_shift_abs(i) / abs(theor_shifts(i))) * 100;
                end
                
                % Calculate refractive index variation (Delta_n) metrics using matching sign conventions
                dn_exp(i) = (-all_pks(best_idx) * n_ave) / nu0;
                err_dn_abs(i) = abs(dn_exp(i) - theor_delta_n(i));
                if theor_delta_n(i) ~= 0
                    err_dn_rel(i) = (err_dn_abs(i) / abs(theor_delta_n(i))) * 100;
                end
            end
        end
    end
    
    % COMPARATIVE TABLE 1: LOCALIZATION (Theoretical vs Experimental positions)
    fprintf('%s\n', repmat('=', 1, 100));
    fprintf('                                    COMPARISON TABLE - LOCALIZATION\n');
    fprintf('%s\n', repmat('=', 1, 100));
    fprintf('%-8s | %-16s | %-16s | %-15s | %-14s | %-18s\n', 'Event', 'Theoretical (m)', ...
        'Experimental (m)', 'Err. Center (m)', 'Err. Width (m)', 'Rel. Error (%)');
    fprintf('%s\n', repmat('-', 1, 100));
    for i = 1:num_events
        str_teo_loc = sprintf('[%.2f; %.2f]', theor_starts(i), theor_ends(i));
        if ~is_detected(i)
            fprintf('Event %-2d | %-16s | %-16s | %-15s | %-14s | %-18s\n', i, ...
                str_teo_loc, 'NOT DETECTED', 'N/A', 'N/A', 'N/A');
        else
            fprintf('Event %-2d | %-16s | %-16s | %-15.4f | %-14.4f | %-18.2f\n', ...
                i, str_teo_loc, str_det_loc{i}, err_center(i), err_width(i), err_width_rel(i));
        end
    end
    fprintf('%s\n\n\n', repmat('=', 1, 100));

    
    % COMPARATIVE TABLE 2: FREQUENCY SHIFT (Theoretical vs Experimental MHz shifts)
    fprintf('%s\n', repmat('=', 1, 77));
    fprintf('                       COMPARISON TABLE - FREQUENCY SHIFT\n');
    fprintf('%s\n', repmat('=', 1, 77));
    fprintf('%-8s | %-13s | %+-12s | %-16s | %-14s\n', 'Event', ...
        'Theor. (MHz)', 'Exper. (MHz)', 'Abs. Error (MHz)', 'Rel. Error (%)');
    fprintf('%s\n', repmat('-', 1, 77));
    for i = 1:num_events
        if ~is_detected(i)
            fprintf('Event %-2d | %+-13.2f | %-12s | %-16s | %-14s\n', i, ...
                theor_shifts(i), 'NOT DETECTED', 'N/A', 'N/A');
        else
            fprintf('Event %-2d | %+-13.2f | %+-12.2f | %-16.2f | %-14.2f\n', ...
                i, theor_shifts(i), shift_MHz_exp(i), err_shift_abs(i), ...
                err_shift_rel(i));
        end
    end
    fprintf('%s\n\n\n', repmat('=', 1, 77));
    
    % COMPARATIVE TABLE 3: DELTA_N (Theoretical vs Calculated refractive index variations)
    fprintf('%s\n', repmat('=', 1, 68));
    fprintf('                     COMPARISON TABLE - DELTA N\n');
    fprintf('%s\n', repmat('=', 1, 68));
    fprintf('%-8s | %-10s | %-12s | %-11s | %-14s\n', 'Event', 'Theor.', ...
        'Exper.', 'Abs. Error', 'Rel. Error (%)');
    fprintf('%s\n', repmat('-', 1, 68));
    for i = 1:num_events
        if ~is_detected(i)
            fprintf('Event %-2d | %-10.2e | %-12s | %-11s | %-14s\n', i, ...
                theor_delta_n(i), 'NOT DETECTED', 'N/A', 'N/A');
        else
            fprintf('Event %-2d | %+-10.2e | %+-12.2e | %-11.2e | %-14.2f\n', ...
                i, theor_delta_n(i), dn_exp(i), err_dn_abs(i), err_dn_rel(i));
        end
    end
    fprintf('%s\n', repmat('=', 1, 68));
end

function [idx_start, idx_end] = find_event_bounds(smooth_freq_shift, idx_peak, threshold)
    % Initialize boundaries at the localized peak index
    idx_start = idx_peak;
    
    % Regress backward down the spatial trace until the frequency shift 
    % drops below 30% of the target noise floor threshold
    while idx_start > 1 && abs(smooth_freq_shift(idx_start)) > (threshold * 0.3)
        idx_start = idx_start - 1;
    end
    
    % Progress forward down the spatial trace until the frequency shift
    % drops below 30% of the target noise floor threshold
    idx_end = idx_peak;
    while idx_end < length(smooth_freq_shift) && abs(smooth_freq_shift(idx_end)) > (threshold * 0.3)
        idx_end = idx_end + 1;
    end
end
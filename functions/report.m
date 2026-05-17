function report(theor_starts, theor_ends, theor_shifts, theor_delta_n, ...
                               all_locs, all_pks, z_valid, smooth_freq_shift, ...
                               threshold, n_ave, nu0)

    num_events = length(theor_starts);
    num_detected = length(all_locs);

    % ========================================================================
    % --- SECTION 7.5: ISOLATED DETECTION REPORT ---
    % ========================================================================
    fprintf('\n=================================================================\n');
    fprintf('                   DETECTION REPORT - PHI-OTDR                      \n');
    fprintf('=================================================================\n');
    fprintf('%-8s | %-16s | %-16s | %-15s\n', 'Event', 'Location (m)', ...
        'Freq. Shift (MHz)', 'Delta_n');
    fprintf('%s\n', repmat('-', 1, 65));
    
    for i = 1:num_detected
        [~, idx_peak] = min(abs(z_valid - all_locs(i)));
        
        idx_start = idx_peak;
        while idx_start > 1 && abs(smooth_freq_shift(idx_start)) > ...
            (threshold * 0.3)
            idx_start = idx_start - 1;
        end
        idx_end = idx_peak;
        while idx_end < length(z_valid) && abs(smooth_freq_shift(idx_end)) > ...
            (threshold * 0.3)
            idx_end = idx_end + 1;
        end
        
        z_start = z_valid(idx_start);
        z_end   = z_valid(idx_end);
        loc_str = sprintf('[%.2f; %.2f]', z_start, z_end);
        
        % De acordo com a tua secção 7.5 original:
        shift_Hz = -all_pks(i); 
        dn_calc_val = (shift_Hz * n_ave) / nu0;
        
        fprintf('Event %-2d | %-16s | %+-16.2f | %+-15.2e\n', ...
            i, loc_str, shift_Hz/1e6, dn_calc_val);
    end
    fprintf('%s\n', repmat('-', 1, 65));

    % ========================================================================
    % --- STATISTICAL PERFORMANCE ANALYSIS (TP, FP, FN) ---
    % ========================================================================
    tolerance = 2.0; % Spatial tolerance margin (meters)
detected_matched = false(1, num_detected); 
    TP = 0; 
    FN = 0; 
    
    fprintf('\n%s\n', repmat('=', 1, 65));
    fprintf('                     EVENT MATCHING ANALYSIS                               \n');
    fprintf('%s\n', repmat('=', 1, 65));
    
    % Step 1: Evaluate True Positives (TP) and False Negatives (FN)
    for i = 1:num_events
        limit_start = theor_starts(i) - tolerance;
        limit_end   = theor_ends(i) + tolerance;
        
        match_idx = find(all_locs >= limit_start & all_locs <= limit_end);
        
        if ~isempty(match_idx)
            TP = TP + 1;
            detected_matched(match_idx) = true; 
            fprintf('TP (Event %d at [%.2f; %.2f] m) :: DETECTED at peak %.2f m\n', ...
                i, theor_starts(i), theor_ends(i), all_locs(match_idx(1)));
        else
            FN = FN + 1;
            fprintf('FN (Event %d at [%.2f; %.2f] m) :: NOT DETECTED by the sensor\n', ...
                i, theor_starts(i), theor_ends(i));
        end
    end
    
    % Step 2: Evaluate False Positives (FP)
    FP = sum(~detected_matched);
    if FP > 0
        fp_indices = find(~detected_matched);
        for k = 1:length(fp_indices)
            fprintf('FP (Noise Artifact) :: DETECTED at peak %.2f m\n', all_locs(fp_indices(k)));
        end
    else
        fprintf('No False Alarms (FP) detected outside the perturbation zones.\n');
    end
    
    % Calculations for percentages and classification metrics
    pct_TP = (TP / num_events) * 100;
    pct_FN = (FN / num_events) * 100;
    if num_detected > 0
        pct_FP = (FP / num_detected) * 100;
    else
        pct_FP = 0;
    end
    
    Sensitivity = (TP / (TP + FN)) * 100; 
    if (TP + FP) > 0
        Precision = (TP / (TP + FP)) * 100;
    else
        Precision = 0;
    end
    
    fprintf('%s\n', repmat('=', 1, 65));
    fprintf('                       SUMMARIZE REPORT                               \n');
    fprintf('%s\n', repmat('=', 1, 65));
    fprintf('Total Real Perturbations:    %d\n', num_events);
    fprintf('Total Detected Peaks:        %d\n', num_detected);
    fprintf('%s\n', repmat('-', 1, 65));
    fprintf('True Positives (TP):         %d  (%.2f%% of real events detected)\n', TP, pct_TP);
    fprintf('False Negatives (FN):        %d  (%.2f%% of real events missed)\n', FN, pct_FN);
    fprintf('False Positives (FP):        %d  (%.2f%% of triggers were false alarms)\n', FP, pct_FP);
    fprintf('%s\n', repmat('-', 1, 65));
    fprintf('System Sensitivity (Recall): %.1f%%\n', Sensitivity);
    fprintf('System Precision:            %.1f%%\n', Precision);
    fprintf('%s\n', repmat('=', 1, 65));

    % ========================================================================
    % --- SECTION 7.6: DIRECT COMPARATIVE TRANSLATION PROCESSING ---
    % ========================================================================
    err_center = NaN(num_events, 1);
    err_length = NaN(num_events, 1);
    str_det_loc = cell(num_events, 1);
    shift_detected_MHz = zeros(num_events, 1);
    err_shift_abs = NaN(num_events, 1);
    err_shift_rel = NaN(num_events, 1);
    dn_calc = zeros(num_events, 1);
    err_dn_abs = NaN(num_events, 1);
    err_dn_rel = NaN(num_events, 1);
    is_detected = false(num_events, 1);
    
    for i = 1:num_events
        t_center = (theor_starts(i) + theor_ends(i)) / 2;
        t_length = theor_ends(i) - theor_starts(i);
        
        if ~isempty(all_locs)
            [dist_min, best_idx] = min(abs(all_locs - t_center));
            
            if dist_min < 4
                is_detected(i) = true;
                peak_loc = all_locs(best_idx);
                [~, idx_peak] = min(abs(z_valid - peak_loc));
                
                idx_start = idx_peak;
                while idx_start > 1 && abs(smooth_freq_shift(idx_start)) > (threshold * 0.3)
                    idx_start = idx_start - 1;
                end
                idx_end = idx_peak;
                while idx_end < length(z_valid) && abs(smooth_freq_shift(idx_end)) > (threshold * 0.3)
                    idx_end = idx_end + 1;
                end
                
                z_start = z_valid(idx_start);
                z_end   = z_valid(idx_end);
                str_det_loc{i} = sprintf('[%.2f; %.2f]', z_start, z_end);
                
                % 1. Spatial Localization Metrics
                exp_center = (z_start + z_end) / 2;
                exp_length = z_end - z_start;
                err_center(i) = abs(exp_center - t_center);
                err_length(i) = exp_length - t_length;
                
                % 2. Frequency Shift Metrics
                shift_Hz_corrected = -all_pks(best_idx);
                shift_detected_MHz(i) = shift_Hz_corrected / 1e6; 
                
                err_shift_abs(i) = abs(shift_detected_MHz(i) - theor_shifts(i));
                if theor_shifts(i) ~= 0
                    err_shift_rel(i) = (err_shift_abs(i) / abs(theor_shifts(i))) * 100;
                else
                    err_shift_rel(i) = 0;
                end
                
                % 3. Refractive Index Variation Metrics (Delta_n)
                dn_calc(i) = (shift_Hz_corrected * n_ave) / nu0;
                err_dn_abs(i) = abs(dn_calc(i) - theor_delta_n(i));
                if theor_delta_n(i) ~= 0
                    err_dn_rel(i) = (err_dn_abs(i) / abs(theor_delta_n(i))) * 100;
                else
                    err_dn_rel(i) = 0;
                end
            end
        end
    end
    
    % --- CROSS-COMPARATIVE TABLES PRINTING ---
    sep_line = repmat('=', 1, 86);
    sub_line = repmat('-', 1, 86);
    
    % COMPARATIVE TABLE 1: LOCALIZATION
    fprintf('\n%s\n                         COMPARISON TABLE - PHI-OTDR LOCALIZATION                        \n%s\n', sep_line, sep_line);
    fprintf('\n                                      Location\n%s\n', sub_line);
    fprintf('%-8s | %-16s | %-16s | %-16s | %-14s\n', 'Event', 'Theor. (m)', ...
        'Exper. (m)', '|Err. Center| (m)', 'Err. Length (m)');
    fprintf('%s\n', sub_line);
    for i = 1:num_events
        str_teo_loc = sprintf('[%.2f; %.2f]', theor_starts(i), theor_ends(i));
        if ~is_detected(i)
            fprintf('Event %-2d | %-16s | %-16s | %-17s | %-14s\n', i, ...
                str_teo_loc, 'NOT DETECTED', 'N/A', 'N/A');
        else
            fprintf('Event %-2d | %-16s | %-16s | %-17.2f | %-+13.2f\n', ...
                i, str_teo_loc, str_det_loc{i}, err_center(i), err_length(i));
        end
    end
    
    % COMPARATIVE TABLE 2: FREQUENCY SHIFT
    fprintf('%s\n\n                                     Frequency Shift\n%s\n', ...
        sep_line, sub_line);
    fprintf('%-8s | %-13s | %+-13s | %-16s | %-14s\n', 'Event', ...
        'Theor. (MHz)', 'Exper. (MHz)', 'Abs. Error (MHz)', 'Rel. Error (%)');
    fprintf('%s\n', sub_line);
    for i = 1:num_events
        if ~is_detected(i)
            fprintf('Event %-2d | %+-13.2f | %-13s | %-16s | %-14s\n', i, ...
                theor_shifts(i), 'NOT DETECTED', 'N/A', 'N/A');
        else
            fprintf('Event %-2d | %+-13.2f | %+-13.2f | %-16.2f | %-14.2f\n', ...
                i, theor_shifts(i), shift_detected_MHz(i), err_shift_abs(i), ...
                err_shift_rel(i));
        end
    end
    
    % COMPARATIVE TABLE 3: DELTA_N
    fprintf('%s\n\n                                          Delta_n\n%s\n', ...
        sep_line, sub_line);
    fprintf('%-8s | %-10s | %-12s | %-11s | %-14s\n', 'Event', 'Theor.', ...
        'Exper.', 'Abs. Error', 'Rel. Error (%)');
    fprintf('%s\n', sub_line);
    for i = 1:num_events
        if ~is_detected(i)
            fprintf('Event %-2d | %-10.2e | %-12s | %-11s | %-14s\n', i, ...
                theor_delta_n(i), 'NOT DETECTED', 'N/A', 'N/A');
        else
            fprintf('Event %-2d | %+-10.2e | %+-12.2e | %-11.2e | %-14.2f\n', ...
                i, theor_delta_n(i), dn_calc(i), err_dn_abs(i), err_dn_rel(i));
        end
    end
    fprintf('%s\n', sep_line);
end
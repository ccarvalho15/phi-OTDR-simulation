function [all_locs, all_pks, threshold] = peak_detection(freq_shift_calib, ...
    smooth_freq_shift, z_valid, conf)
    %% ========================================================================
    % 7. PEAK DETECTION WITH NOISE SUPPRESSION
    % ========================================================================
    %   IN : freq_shift, freq_shift_calib, z_valid, M, dz, pert_length, std_mult
    %   OUT: locs_pos, pks_pos — positions and amplitudes of positive peaks
    %        locs_neg, pks_neg — positions and amplitudes of negative peaks
    %        all_locs, all_pks — merged, sorted event list
    %        smooth_freq_shift — low-pass smoothed shift profile
    %        threshold         — adaptive noise floor
    
    std_mult = conf.std_mult;
    pert_length = conf.pert_length;
    dz = conf.dz;
    
    % -------------------------------------------------------------------
    % 7.2 ADAPTIVE NOISE THRESHOLD
    % -------------------------------------------------------------------
    %   The standard deviation of the calibration trace reflects the noise
    %   floor; events must exceed std_mult*sigma(freq_shift_calib) to be
    %   declared.
    threshold = std_mult * std(freq_shift_calib);


    % -------------------------------------------------------------------
    % 7.3 SEPARATED POSITVE/NEGATIVE PEAK SEARCH
    % -------------------------------------------------------------------
    %   Events with opposite-sign delta_n produce opposite-sign shifts, so
    %   positive and negative peaks are detected independently.
    min_peak_dist = round((pert_length*1.5) / dz); % Minimum inter-peak grid spacing

    % Common findpeaks settings:
    %   MinPeakHeight :: rejects all variations below the threshold floor
    %   MinPeakDistance :: ensures detected peaks are separated by at least one
    %                      event width
    %   MinPeakProminence :: measures peak height relative to the surrounding
    %                        baseline, filtering out ripple artifacts located
    %                        on long signal slopes
    peak_opts = { ...
        'MinPeakHeight',      threshold, ...
        'MinPeakDistance',    min_peak_dist, ...
        'MinPeakProminence',  threshold * 0.9};

    [pks_pos, idx_pos] = findpeaks( smooth_freq_shift, peak_opts{:});
    [pks_neg, idx_neg] = findpeaks(-smooth_freq_shift, peak_opts{:});
    pks_neg = -pks_neg; % Restore true (negative) amplitudes
    locs_pos = z_valid(idx_pos);
    locs_neg = z_valid(idx_neg);

    % -------------------------------------------------------------------
    % 7.4 MERGE AND STOR ALL DETECTED EVENTS BY POSITION
    % -------------------------------------------------------------------
    all_locs = [locs_pos, locs_neg];
    all_pks  = [pks_pos,  pks_neg];
    [all_locs, sort_idx] = sort(all_locs);
    all_pks = all_pks(sort_idx);
end
function [freq_shift, freq_shift_calib, smooth_freq_shift, corr_map, lags_freq] = ...
    spectral_shift_estimation(conf, E_ref, E_sig, E_ref_raw_all)

    %% ========================================================================
    % 6. SPECTRAL SHIFT ESTIMATION VIA CROSS-CORRELATION
    % ========================================================================
    %   IN : E_ref, E_sig, E_ref_raw_all, delta_f, Nf, Nz, M, SNR_dB
    %   OUT: freq_shift(1×Nz-M+1) :: per-position spectral shift estimate
    %        freq_shift_calib(1×Nz-M+1) :: calibration trace (noise vs noise)
    %        corr_map(2Nf-1 × Nz-M+1) :: full correlation surface for 3-D plot
    %        lags_freq :: frequency lag axis (Hz)
    %
    %   At each position k, the intensity spectra |E|^2 of the reference and
    %   perturbed states are cross-correlated. The lag that maximises the
    %   normalised cross-correlation coefficient is the local spectral shift
    %   Delta_nu, which encodes the environmental change at that point.
    %
    %   Zero-mean spectra are used to eliminate DC bias; 'coeff' normalises the
    %   output so the auto-correlation peak at zero lag equals 1.
    % 
    %   Following the methodology described in the paper, we perform a local 
    %   cross-correlation between the reference and perturbed Rayleigh 
    %   backscatter Intensity spectra to quantify the environmental impact.
    
    Nf = conf.Nf;
    Nz = conf.Nz;
    delta_f = conf.delta_f;
    M = conf.M;
    SNR_dB = conf.SNR_dB;
    
    lags_freq = (-(Nf-1):(Nf-1)) * delta_f; % Lag axis in Hz
    % Pre-allocate matrix for 3D correlation visualization
    corr_map = zeros(length(lags_freq), Nz-M+1);
    % Vector to store the estimated frequency shift per position
    freq_shift = zeros(1, Nz-M+1);               
    freq_shift_calib = zeros(1, Nz - M +1);
    
    for k = 1:Nz-M+1
        % Intensity spectra (power spectral density approximation)
        I_ref = abs(E_ref(:,k)).^2;
        I_sig = abs(E_sig(:,k)).^2;
        
        % Normalised cross-correlation of zero-mean spectra
        [cv, lags] = xcorr(I_sig - mean(I_sig), I_ref - mean(I_ref), 'coeff');
        corr_map(:, k) = cv;
        
        % Peak lag -> spectral shift at position k
        [~, max_idx] = max(cv);
        freq_shift(k) = lags(max_idx) * delta_f; 
    
        % Calibration path: correlate two independent realisations of the
        % reference state to estimate the noise-induced shift variance
    
        I_calib = abs(awgn(E_ref_raw_all(:,k), SNR_dB, 'measured')).^2;
        [cv_c, lags_c] = xcorr(I_calib - mean(I_calib), I_ref - mean(I_ref), 'coeff');
        [~, max_idx_c] = max(cv_c);
        freq_shift_calib(k) = lags_c(max_idx_c) * delta_f;
    end

    % -------------------------------------------------------------------
    % 7.1 SPATIAL SMOOTHING
    % -------------------------------------------------------------------
    %   A moving average with window M suppresses high-frequency speckle noise
    %   without broadening real events (whose spatial extent >> M·dz).
    %
    %   'movmean' computes a moving average across adjacent spatial points. 
    %   This filters out Isolated high-frequency noise spikes without 
    %   distorting broad real events.
    %   Why avoid the raw signal?
    %       Running findpeaks on the raw trace generates hundreds of false 
    %       positives from the random +/-10 MHz noise spikes. 
    
    smooth_window = M;
    smooth_freq_shift = movmean(freq_shift, smooth_window);

end
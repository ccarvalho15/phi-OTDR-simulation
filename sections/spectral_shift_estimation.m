function [freq_shift, freq_shift_calib, smooth_freq_shift, corr_map, lags_freq] = ...
    spectral_shift_estimation(conf, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_calib_all)

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
    z = conf.z;

    z_valid = z(1 : (Nz - M + 1));
    
    lags_freq = (-(Nf-1):(Nf-1)) * delta_f; % Lag axis in Hz
    % Pre-allocate matrix for 3D correlation visualization
    corr_map = zeros(length(lags_freq), Nz-M+1);
    corr_map_id = zeros(length(lags_freq), Nz-M+1);
    % Vector to store the estimated frequency shift per position
    freq_shift = zeros(1, Nz-M+1);      
    freq_shift_id = zeros(1, Nz-M+1);
    freq_shift_calib = zeros(1, Nz - M +1);
    
    for k = 1:Nz-M+1
        % Intensity spectra (power spectral density approximation)
        I_ref = abs(E_ref(:,k)).^2;
        I_sig = abs(E_sig(:,k)).^2;
        
        % Normalised cross-correlation of zero-mean spectra
        % [cv, lags] = xcorr(I_sig - mean(I_sig), I_ref - mean(I_ref), 'coeff');
        [cv, lags] = xcorr(I_ref - mean(I_ref), I_sig - mean(I_sig), 'coeff');
        corr_map(:, k) = cv;
        
        % Peak lag -> spectral shift at position k
        [~, max_idx] = max(cv);
        freq_shift(k) = lags(max_idx) * delta_f; 

        % IDEAL
        I_ref_id = abs(E_ref_id(:, k)).^2;
        I_sig_id = abs(E_sig_id(:, k)).^2;
        
        [cv_id, ~] = xcorr(I_ref_id - mean(I_ref_id), I_sig_id - mean(I_sig_id), 'coeff');
        corr_map_id(:, k) = cv_id;
        
        [~, max_idx_id] = max(cv_id);
        freq_shift_id(k) = lags_freq(max_idx_id);
    
        % Calibration path: correlate two independent realisations of the
        % reference state to estimate the noise-induced shift variance
    
        I_calib = abs(awgn(E_ref_raw_calib_all(:,k), SNR_dB, 'measured')).^2;
        [cv_c, lags_c] = xcorr(I_calib - mean(I_calib), I_ref - mean(I_ref), 'coeff');
        [~, max_idx_c] = max(cv_c);
        freq_shift_calib(k) = lags_c(max_idx_c) * delta_f;
    end

    % %% ========================================================================
    % % FIGURA 5: ESPECTROGRAMA / MAPA 3D IDEAL (SEM RUÍDO)
    % % ========================================================================
    % fig5 = figure(5); set(fig5, 'Name', 'Cross-Correlation Ideal');
    % [Z_mesh, F_mesh] = meshgrid(z_valid(1:Nz-M+1), lags_freq / 1e6);
    % surf(Z_mesh, F_mesh, corr_map_id, 'EdgeColor', 'none');
    % view(35, 45); colormap('jet'); cb = colorbar; cb.FontSize = 12;
    % set(gca, 'FontSize', 15, 'LineWidth', 1.5);
    % xlabel('Distance (m)', 'FontSize', 16); 
    % ylabel('Frequency Lag (MHz)', 'FontSize', 16);
    % zlabel('Correlation', 'FontSize', 16);
    % title('Ideal Spectral Cross-Correlation Map (Noise-Free Baseline)', ...
    %     'FontSize', 18, 'FontWeight', 'bold');
    % rotate3d on;
    % xlim([200 240]); ylim([-1100 1100]); zlim([-0.5 1]);
    % set(gca, 'Position', [0.12, 0.15, 0.65, 0.72]); % Proporção padrão
    % cb.Position = [0.84, 0.15, 0.025, 0.72];
    % 
    % %% ========================================================================
    % % FIGURA 6: ESPECTROGRAMA / MAPA 3D REAL (COM RUÍDO)
    % % ========================================================================
    % fig6 = figure(6); set(fig6, 'Name', 'Cross-Correlation Real');
    % [Z_mesh, F_mesh] = meshgrid(z_valid(1:Nz-M+1), lags_freq / 1e6);
    % surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none');
    % view(35, 45); colormap('jet'); cb = colorbar; cb.FontSize = 12;
    % set(gca, 'FontSize', 16, 'LineWidth', 1.5);
    % xlabel('Distance (m)', 'FontSize', 16); 
    % ylabel('Frequency Lag (MHz)', 'FontSize', 16);
    % zlabel('Correlation', 'FontSize', 16);
    % title(sprintf('Cross-Correlation Map under Noisy Conditions (SNR = %d dB)', ...
    %     SNR_dB), 'FontSize', 18, 'FontWeight', 'bold');
    % rotate3d on;
    % xlim([200 240]); ylim([-1100 1100]); zlim([-0.5 1]);
    % set(gca, 'Position', [0.12, 0.15, 0.65, 0.72]); % Proporção padrão
    % cb.Position = [0.84, 0.15, 0.025, 0.72];

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

    % smooth_small = movmean(freq_shift, max(1, round(M/2)));
    % smooth_large = movmean(freq_shift, 2*M);

    % % Criação da figura
    % fig7 = figure(7);
    % set(fig7, 'Name', 'Raw vs Smooth Frequency Shift');
    % 
    % % 1. Linha do sinal NORMAL (bruto, em vermelho/laranja com alguma transparência)
    % plot(z_valid, freq_shift / 1e6, 'r-', ...
    %     'LineWidth', 1.5, 'DisplayName', 'Raw Trace', 'Color', [0.4 0.4 0.4]);
    % hold on;
    % % 2. Linha do sinal SMOOTH (suavizado, em azul forte e mais espesso)
    % plot(z_valid, smooth_freq_shift / 1e6, 'b-', ...
    %     'LineWidth', 1.5, 'DisplayName', 'Smoothed Trace', 'Color', [0.2 0.6 1.0]);
    % hold off;
    % 
    % % Ajustes de eixos, grelha e legenda
    % grid on;
    % set(gca, 'FontSize', 15, 'LineWidth', 1.3);
    % xlabel('Distance (m)', 'FontSize', 16);
    % ylabel('Frequency Shift (MHz)', 'FontSize', 16);
    % title('Distributed Spectral Shift Profile (Raw vs. Smooth)', ...
    %     'FontSize', 18, 'FontWeight', 'bold');
    % legend('Location', 'northeast', 'FontSize', 14);
    % xlim([205, 240]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Criação da Figura 8
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fig8 = figure(8);
    % set(fig8, 'Name', 'Smoothing Window Trade-off Analysis');
    % 
    % % 1. Sinal Bruto (Raw Trace)
    % plot(z_valid, freq_shift / 1e6, '-', ...
    %     'LineWidth', 1.0, 'DisplayName', 'Raw Trace', 'Color', [0.7 0.7 0.7]);
    % hold on;
    % 
    % % 2. Suavização Insuficiente (M/2) - Deixa passar speckle
    % plot(z_valid, smooth_small / 1e6, '--', ...
    %     'LineWidth', 1.4, 'DisplayName', sprintf('Under-smoothed (M = %d)', max(1, round(M/2))), ...
    %     'Color', [0.85 0.32 0.09]); % Laranja
    % 
    % % 3. Suavização Ideal (M) - Equilibra speckle e resolução espacial
    % plot(z_valid, smooth_freq_shift / 1e6, '-', ...
    %     'LineWidth', 2.0, 'DisplayName', sprintf('Optimal (M = %d)', smooth_window), ...
    %     'Color', [0.00 0.45 0.74]); % Azul principal
    % 
    % % 4. Suavização Excessiva (2M) - Alarga as transições dos eventos
    % plot(z_valid, smooth_large / 1e6, '-.', ...
    %     'LineWidth', 1.6, 'DisplayName', sprintf('Over-smoothed (M = %d)', 2*M), ...
    %     'Color', [0.47 0.67 0.19]); % Verde
    % hold off;
    % 
    % % Ajustes Visuais
    % grid on;
    % set(gca, 'FontSize', 15, 'LineWidth', 1.3);
    % xlabel('Distance (m)', 'FontSize', 16);
    % ylabel('Frequency Shift (MHz)', 'FontSize', 16);
    % title('Trade-off of Moving-Average Window Size on Spectral Shift Profile', ...
    %     'FontSize', 17, 'FontWeight', 'bold');
    % legend('Location', 'northeast', 'FontSize', 13);
    % xlim([205, 240]);

end
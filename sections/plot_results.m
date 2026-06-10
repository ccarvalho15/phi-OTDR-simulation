function plot_results(conf, z_valid, freq_shift, smooth_freq_shift, ...
    corr_map, lags_freq, all_locs, all_pks, threshold)

    z     = conf.z;
    Nz    = conf.Nz;
    M     = conf.M;
    L     = conf.L;

    locs_pos = all_locs(all_pks > 0);
    pks_pos  = all_pks(all_pks > 0);
    locs_neg = all_locs(all_pks < 0);
    pks_neg  = all_pks(all_pks < 0);
    
    first_event = conf.sensing_zone;
    last_event  = conf.sensing_zone + (conf.num_events-1) * ...
        (conf.pert_length + conf.spacing) + conf.pert_length;
    
    % --- FIGURE 1: 1D Frequency Shift Profile ---
    % This plot represents the demodulated sensing signal.
    % The peaks should align with the 'start_idx' and 'end_idx' defined in Section 3.
    
    figure(2)
    set(gcf,  'Name', 'Freq Shift Profile');
    plot(z(1:Nz-M+1), freq_shift / 1e6, 'LineWidth', 1.5)
    grid on; ylabel('Frequency Shift (MHz)'); xlabel('Distance (m)');
    title('Detected Frequency Shift along the Fiber');
    
    
    %%
    % --- FIGURE 2: 3D Correlation Surface ---
    % Provides a global view of the cross-correlation peaks.
    % High correlation (close to 1.0) indicates high similarity between
    % the reference and signal Rayleigh spectra at the shifted frequency.
    
    % figure(3)
    % set(gcf,  'Name', '3D CC (reduced)');
    % [Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
    % surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none')
    % view(35, 45); colormap('jet'); colorbar;
    % xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)'); zlabel('Correlation');
    % title('3D Cross-Correlation Map');
    % rotate3d on;
    % xlim([max(210, first_event - 10) min(L, last_event + 10)]);  % Fiber length
    % ylim([-250 250]); % Frequency shift window
    % zlim([-0.5 1]); % Correlation magnitude scale
    
    
    %%
    % --- FIGURE 3: 2D Contour Map with Peak Tracking ---
    % A "top-down" view that combines the correlation energy with the
    % mathematical peak detection (white line).
    % This visualization is excellent for assessing the Signal-to-Noise Ratio (SNR).
    
    figure(4)
    set(gcf, 'Name', '2D C. Peak Tracking');
    [Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
    contourf(Z_mesh, F_mesh, corr_map, 20, 'LineColor', 'none');
    colormap('jet'); colorbar;
    hold on;
    % The white line traces the maximum correlation lag, verifying the algorithm's accuracy
    plot(z(1:Nz-M+1), freq_shift / 1e6, 'w', 'LineWidth', 1.5);
    hold off;
    xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)');
    title('2D Correlation Map (Top View with Peak Trace)');
    % Focus the view on the perturbed regions for better detail
    xlim([max(0, first_event - 50) min(L, last_event + 50)]);
    ylim([-250 250]);
    
    
    %%
    % --- FIGURE 4: Summary Overview ---
    % Combined plot for comparative analysis of spatial and spectral data.
    figure(5)
    set(gcf, 'Name', 'Overview');
    subplot(2,1,1)
    plot(z(1:Nz-M+1), freq_shift / 1e6, 'LineWidth', 1.5)
    grid on; ylabel('Frequency Shift (MHz)'); xlabel('Distance (m)');
    title('Detected Frequency Shift along the Fiber');
    
    subplot(2,1,2)
    [Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
    surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none')
    view(35, 45); colormap('jet'); colorbar;
    xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)'); zlabel('Correlation');
    title('3D Cross-Correlation Map');
    rotate3d on;
    
    xlim([0 L]); % Full fiber length
    ylim([-250 250]); % Frequency shift window
    zlim([-0.5 1]); % Correlation magnitude scale
    
    
    %%
    % --- FIGURE 5: ATTENUATION IMPACT ANALYSIS -----
    
    % P_ideal_mean = mean(abs(E_ref_id).^2, 1);
    % P_noisy_mean = mean(abs(E_ref).^2,    1);
    %
    % % Normalizing electric field
    % % Calculates the scaling factor to map dimensionless simulated fields to
    % % physical mW. Using mean of the first point to stabilize against coherent
    % % fading (speckle).
    % scale_factor = P_input_mW / max(P_ideal_mean);
    %
    % % Select the first frequency from the sweep for visualization
    % P_ref_ideal = P_ideal_mean * scale_factor;
    % P_ref_loss  = P_noisy_mean * scale_factor;
    %
    % % Comparing the backscattered intensity with and without fiber loss
    % figure(6)
    % set(gcf, 'Name', 'Attenuation Impact');
    %
    % % Subplot 2: Logarithmic Scale (OTDR Trace)
    % plot(z_valid, 10*log10(P_ref_ideal + eps), 'b', 'DisplayName', ...
    %     'Ideal (No Loss)');
    % hold on;
    % plot(z_valid, 10*log10(P_ref_loss + eps), 'r', 'DisplayName', ...
    %     ['Fiber Loss (', num2str(attenuation), ' dB/km) + Noise']);
    % hold off;
    % title('Backscattered Intensity (Logarithmic Scale)');
    % ylabel('Power (dBm)'); xlabel('Distance (m)'); grid on;  legend('Location', 'southeast');
    
    %%
    % --- FIGURE 7: Peak Detection Diagnostics ---
    % Visualizes raw versus smoothed frequency shifts alongside thresholds and
    % positive/negative peak markers.
    figure(7);
    
    plot(z_valid, freq_shift/1e6, 'Color', [0.4 0.4 0.4], 'DisplayName', 'Raw shift');
    hold on;
    plot(z_valid, smooth_freq_shift/1e6, 'Color', [0.2 0.6 1.0], 'LineWidth', 1.5, ...
        'DisplayName', 'Smoothed shift');
    scatter(locs_pos, pks_pos/1e6, 80, 'g^', 'filled', 'DisplayName', 'Detected (+)');
    scatter(locs_neg, pks_neg/1e6, 80, 'rv', 'filled', 'DisplayName', 'Detected (−)');
    yline( threshold/1e6, 'g--', 'LineWidth', 1, 'DisplayName', '+Threshold');
    yline(-threshold/1e6, 'r--', 'LineWidth', 1, 'DisplayName', '−Threshold');
    hold off;
    xlabel('Distance (m)'); ylabel('Frequency Shift (MHz)');
    title('Peak Detection — φ-OTDR');
    legend('Location', 'northwest'); grid on;
    xlim([0 240])


end
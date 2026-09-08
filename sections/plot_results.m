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
    
    %%
    % --- FIGURE 1: 1D Frequency Shift Profile ---
    % This plot represents the demodulated sensing signal.
    % The peaks should align with the 'start_idx' and 'end_idx' defined in Section 3.
    
    fig2 = figure(2);
    set(fig2, 'Name', '1D Frequency Shift Profile');
    plot(z(1:Nz-M+1), freq_shift / 1e6, 'LineWidth', 1.5);
    grid on; 
    set(gca, 'FontSize', 15, 'LineWidth', 1.5);
    ylabel('Frequency Shift (MHz)', 'FontSize', 16); 
    xlabel('Distance (m)', 'FontSize', 16);
    title('Distributed Frequency Shift Trace Along Sensing Fiber', ...
        'FontSize', 18, 'FontWeight', 'bold');
    xlim([0 245])
    
    % filename2 = fullfile(output_dir, sprintf('Frequency_Shift_%s.png', timestamp));
    % exportgraphics(fig2, filename2, 'Resolution', 300);
    
    %%
    % --- FIGURE 2: 3D Correlation Surface ---
    % Provides a global view of the cross-correlation peaks.
    % High correlation (close to 1.0) indicates high similarity between
    % the reference and signal Rayleigh spectra at the shifted frequency.
    
    views = {
    [35, 45],   'Perspective View',    '35-45';
    [135, 45],  'Rear View',           '135-45';
    [0, 90],    'Top-Down / 2D View',  '0-90';
    [90, 40],   'Lateral View',        '90-40'
};

    [Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);

    for k = 1:size(views, 1)
        fig_view = figure(10 + k);
        set(fig_view, 'Name', sprintf('3D CC View %s', views{k,3}));
        
        surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none');
        colormap('jet'); 
        
        set(gca, 'FontSize', 15, 'LineWidth', 1.5);
        xlabel('Distance (m)', 'FontSize', 16); 
        ylabel('Frequency Lag (MHz)', 'FontSize', 16); 
        zlabel('Correlation', 'FontSize', 16);
        title(sprintf('Cross-Correlation Surface Map — %s', views{k,2}), ...
          'FontSize', 18, 'FontWeight', 'bold');
        
        xlim([0 240]); 
        ylim([-1100 1100]); 
        zlim([-0.5 1]);
        
        view(views{k,1}(1), views{k,1}(2));
        rotate3d on;

        cb = colorbar;
        cb.FontSize = 12;

        % Ajuste dedicado para a vista 90-40 (evita sobreposição dos rótulos de Distance)

        set(gca, 'Position', [0.12, 0.15, 0.65, 0.72]); % Proporção padrão
        cb.Position = [0.84, 0.15, 0.025, 0.72];
    end
        
    
    %%
    % --- FIGURE 3: 2D Contour Map with Peak Tracking ---
    % A "top-down" view that combines the correlation energy with the
    % mathematical peak detection (white line).
    % This visualization is excellent for assessing the Signal-to-Noise Ratio (SNR).
    
    % figure(4)
    % set(gcf, 'Name', '2D C. Peak Tracking');
    % [Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
    % contourf(Z_mesh, F_mesh, corr_map, 20, 'LineColor', 'none');
    % colormap('jet'); colorbar;
    % hold on;
    % % The white line traces the maximum correlation lag, verifying the algorithm's accuracy
    % plot(z(1:Nz-M+1), freq_shift / 1e6, 'w', 'LineWidth', 1.5);
    % hold off;
    % xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)');
    % title('2D Correlation Map (Top View with Peak Trace)');
    % % Focus the view on the perturbed regions for better detail
    % xlim([0 240])
    % % xlim([max(0, first_event - 10) min(L, last_event + 10)]);
    % ylim([-1000 1000]);
    
    %%
    % --- FIGURE 4: Summary Overview ---
    % Combined plot for comparative analysis of spatial and spectral data.
    % figure(5)
    % set(gcf, 'Name', 'Overview');
    % subplot(2,1,1)
    % plot(z(1:Nz-M+1), freq_shift / 1e6, 'LineWidth', 1.5)
    % grid on; ylabel('Frequency Shift (MHz)'); xlabel('Distance (m)');
    % title('Detected Frequency Shift along the Fiber');
    % 
    % subplot(2,1,2)
    % [Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
    % surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none')
    % view(35, 45); colormap('jet'); colorbar;
    % xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)'); zlabel('Correlation');
    % title('3D Cross-Correlation Map');
    % rotate3d on;
    % 
    % xlim([0 L]); % Full fiber length
    % ylim([-1200 1200]); % Frequency shift window
    % zlim([-0.5 1]); % Correlation magnitude scale
    
    
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
    
    fig4 = figure(4);
    set(fig4, 'Name', 'Peak Detection Diagnostics');
    
    plot(z_valid, freq_shift/1e6, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, 'DisplayName', 'Raw trace');
    hold on;
    plot(z_valid, smooth_freq_shift/1e6, 'Color', [0.2 0.6 1.0], 'LineWidth', 2.5, 'DisplayName', 'Smoothed trace');
    scatter(locs_pos, pks_pos/1e6, 100, 'g^', 'filled', 'DisplayName', 'Detected (+)');
    scatter(locs_neg, pks_neg/1e6, 100, 'rv', 'filled', 'DisplayName', 'Detected (−)');
    yline( threshold/1e6, 'g--', 'LineWidth', 1.8, 'DisplayName', '+Threshold');
    yline(-threshold/1e6, 'r--', 'LineWidth', 1.8, 'DisplayName', '−Threshold');

    offset_y = max(smooth_freq_shift/1e6) * 0.05; 
    for p = 1:length(locs_pos)
        text(locs_pos(p), (pks_pos(p)/1e6) + offset_y, sprintf('%+.1f MHz', pks_pos(p)/1e6), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.1 0.5 0.1]);
    end
    for n = 1:length(locs_neg)
        text(locs_neg(n), (pks_neg(n)/1e6) - offset_y, sprintf('%+.1f MHz', pks_neg(n)/1e6), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.7 0.1 0.1]);
    end
    hold off;
    
    set(gca, 'FontSize', 15, 'LineWidth', 1.5);
    xlabel('Distance (m)', 'FontSize', 16); 
    ylabel('Frequency Shift (MHz)', 'FontSize', 16);
    title('Distributed Spectral Shift Event Identification via Threshold Detection', ...
        'FontSize', 18, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 14); 
    grid on; 
    xlim([205 245]); ylim([-1100 1100])


end
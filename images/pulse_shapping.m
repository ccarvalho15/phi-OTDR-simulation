clc, clear, close
% 1. SYSTEM CONFIGURATION & WAVEGUIDE PROPERTIES

c = 3e8;                            % Speed of light in vacuum (m/s)
lambda0 = 1550e-9;                  % Operating wavelength (m)
nu0 = c/lambda0;                    % Central optical frequency (Hz)
L = 240;                            % Total fiber length used on the lab (m) 
n_ave = 1.456;                      % Average refractive index of the silica fiber
dz = 0.05;                          % Spatial sampling interval (m)
z = 0:dz:L-dz;                      % Distance vector along the fiber 1D model
Nz = length(z);                     % Total number os spatial sampling points
pulse_width = 10e-9;                % Temporal pulse width (10 ns)
d = c * pulse_width / (2 * n_ave);  % Spatial resolution (~1 m)
M = round(d / dz);                  % Number of scattering segments 
                                    % (inhomogeneities) within one pulse

% 1. Ideal Rectangular
rect_pulse = ones(1, M);

% 2. RC Filter Model (fc = 100 MHz)
fc_comp = 100e6; 
RC_comp = 1 / (2 * pi * fc_comp);
t_rc = linspace(0, 2 * pulse_width, 2 * M);
window_rc_full = zeros(1, 2 * M);
for k = 1:2*M
    tk = t_rc(k);
    if tk <= pulse_width
        window_rc_full(k) = 1 - exp(-tk / RC_comp);
    else
        window_rc_full(k) = (1 - exp(-pulse_width / RC_comp)) * exp(-(tk - pulse_width) / RC_comp);
    end
end
window_rc = window_rc_full(1:M) / max(window_rc_full);

% 3. Super-Gaussian (Order N=3)
t_sg = linspace(-pulse_width, pulse_width, M);
window_sg = sgauss(t_sg, pulse_width, 1, 0, 3);
window_sg = window_sg / max(window_sg);

figure('Name', 'Pulse Shape Comparison Analysis');
t_plot_ns = linspace(0, pulse_width * 1e9, M);


plot(t_plot_ns, window_rc, 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'DisplayName', 'RC Filter (Realistic EOM)');
hold on;
plot(t_plot_ns, window_sg, 'Color', [0 0.45 0.74], 'LineWidth', 2, 'DisplayName', 'Super-Gaussian (N=3)');
plot(t_plot_ns, rect_pulse, '--y', 'LineWidth', 1.2, 'DisplayName', 'Ideal Rectangular');


grid on;
xlabel('Time (ns)');
ylabel('Normalized Amplitude');
title('Comparison of Optical Pulse Shapes');
legend('Location', 'southoutside', 'Orientation', 'horizontal');
ylim([0 1.2]);
xlim([0 pulse_width*1e9]);
hold off;

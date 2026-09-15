clc, clear, close

%%
addpath('sections\')
addpath('functions\')
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
M = 250;                            % Number of scattering segments
                                    % (inhomogeneities) within one pulse
          
fig1 = figure(1);
set(fig1, 'Name', 'Pulse Shape Comparison Analysis');
hold on;

% 1. Ideal Rectangular
rect_pulse = ones(1, M);
t_plot_ns = linspace(0, pulse_width * 1e9, M);
t_rect_ns = [-5, 0, 0, pulse_width*1e9, pulse_width*1e9, 20];
v_rect    = [ 0, 0, 1,               1,               0,  0];
plot(t_rect_ns, v_rect, 'r--', 'LineWidth', 2, 'DisplayName', 'Ideal Rectangular Pulse');


% % 2. RC Filter Model (fc = 100 MHz)
% fc = [50e6, 100e6, 200e6]; 
% for fc_comp = fc
%     RC_comp = 1 / (2 * pi * fc_comp);
%     t_rc = linspace(0, 2 * pulse_width, 2 * M);
%     window_rc_full = zeros(1, 2 * M);
%     for k = 1:2*M
%         tk = t_rc(k);
%         if tk <= pulse_width
%             window_rc_full(k) = 1 - exp(-tk / RC_comp);
%         else
% 
%             window_rc_full(k) = (1 - exp(-pulse_width / RC_comp)) * exp(-(tk - pulse_width) / RC_comp);
%         end
%     end
%     window_rc = window_rc_full / max(window_rc_full);
%     plot(t_rc * 1e9, window_rc, 'LineWidth', 2, ...
%         'DisplayName', sprintf('Filtro RC (f_c = %d MHz)', fc_comp/1e6));
% end


% 3. Super-Gaussian (Order N)
pulse_width_ns = pulse_width * 1e9; % 10 ns
t_sg_ns = linspace(-pulse_width_ns, 2*pulse_width_ns, 2*M); % -10 a 20 ns
t_center_ns = pulse_width_ns / 2; % 5 ns

order_N = [1, 3, 5];

for N = order_N
    window_sg = sgauss(t_sg_ns - t_center_ns, pulse_width_ns, 1, 0, N);
    window_sg = window_sg / max(window_sg);
    plot(t_sg_ns, window_sg, 'LineWidth', 2, 'DisplayName', ['Super-Gaussian (p=' num2str(N) ')']);
end

grid on;
xlim([-11, 21]); 
% xlim([-7, 22])
ylim([0, 1.2]);
set(gca, 'FontSize', 20, 'LineWidth', 1.5);
    
xlabel('Relative Time (ns)', 'FontSize', 20); 
ylabel('Normalized Amplitude', 'FontSize', 20);
%title('Comparison of Optical Pulse Shapes', 'FontSize', 25, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 15);
hold off;

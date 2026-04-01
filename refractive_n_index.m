clear; clc; close all

%% ---------------- PARAMETERS ----------------
c = 3e8;                    % Speed of light (m/s)
lambda0 = 1550e-9;          % Central wavelength (m)
nu0 = c/lambda0;            % Central optical frequency (Hz)
n_ave = 1.456;              % Average refractive index of the silica fiber

L = 3000;                   % Total fiber length (m)
dz = 0.1;                   % Spatial sampling interval (m)
z = 0:dz:L-dz;              % Distance vector
Nz = length(z);

pulse_width = 10e-9;        % Pulse width (10 ns)
d = c*pulse_width/(2*n_ave);% Spatial resolution (~1 m)
M = round(d/dz);            % Number of points within one pulse width
alpha = 0;                  % Fiber attenuation (neglected for this simulation)

%% ------------ RANDOM REFRACTIVE INDEX PROFILE -----------
% Simulating the inhomogeneous nature of the fiber core (Rayleigh centers)
sigma_n = 2e-6;                
delta_n = sigma_n*randn(1,Nz);
n = n_ave + delta_n;           % Base refractive index profile

%% ------------ FRESNEL REFLECTION COEFFICIENT ------------
% Calculating local reflection coefficients along the fiber
r = zeros(1,Nz);
for i = 1:Nz-1
    r(i) = (n(i)-n(i+1))/(n(i)+n(i+1));
end

%% ------------- DIRECT REFRACTIVE INDEX PERTURBATION ------
% Creating a manual refraxtive index change vector
delta_n_pert = zeros(1, Nz);

% Apply a refractive index increase of 1e-4 between 1500 m and 1510 m
delta_n_pert(round(1500/dz):round(1510/dz)) = 1e-4;

% Perturbed profile (original + manual change)
n_pert = n + delta_n_pert;

%% SIMULATION OF STATIC MEASUREMENTS
freq_range = 1000e6;        % Frequency scanning range (1000 MHz)
delta_f = 10e6;             % Frequency step (10 MHz)
% Absolute optical frequency vector centered at nu0
f = nu0 + (-freq_range/2 : delta_f : freq_range/2);  
Nf = length(f);

% Pre-allocating matrices for reference and signal fields (Frequency x Distance)
E_ref = zeros(Nf, Nz-M);
E_sig = zeros(Nf, Nz-M);

%% ------------ FREQUENCY SCANNING LOOP ------------
% Simulating the backscattered field for each frequency in the scan
fprintf('Processing frequenct scan...');
for f_idx = 1:Nf
    current_nu = f(f_idx);
    
    % Reference propagation phase
    beta_ref = 2 * pi * n * current_nu /c;
    phi_ref = cumsum(beta_ref)*dz;      % Integrated phase

    % Signal propagation phase (fiber perturbed by delta_n)
    beta_sig = 2 * pi * n_pert * current_nu / c;
    phi_sig = cumsum(beta_sig)*dz;

    for k = 1:Nz-M
        idx = k:(k+M-1);
        % Vector summation of backscattered fields within the pulse volume
        E_ref(f_idx, k) = sum(r(idx) .* exp(1j * 2 * phi_ref(idx)));
        E_sig(f_idx, k) = sum(r(idx) .* exp(1j * 2 * phi_sig(idx)));
    end
end
fprintf('Done!\n');

%% ------------ CROSS-CORRELATION ------------
% Calculating frequency shift via local cross-correlation of Rayleigh spectra

lags_freq = (-(Nf-1):(Nf-1)) * delta_f; % Frequency lag axis
corr_map = zeros(length(lags_freq), Nz-M); % Matrix for 3D visualization
freq_shift = zeros(1, Nz-M); 

for k = 1:Nz-M
    ref = abs(E_ref(:,k)).^2;
    sig = abs(E_sig(:,k)).^2;
    
    % Cross-correlation of zero-mean intensity spectra
    [cv, lags] = xcorr(sig - mean(sig), ref - mean(ref), 'coeff');
    corr_map(:, k) = cv;
    
    % Finding the peak index to determine the frequency shift
    [~, max_idx] = max(cv);
    freq_shift(k) = lags(max_idx) * delta_f; 
end

%% ------------ 3D VISUALIZATION ------------
figure(1)
plot(z(1:Nz-M), freq_shift / 1e6, 'LineWidth', 1.5)
grid on; ylabel('Frequency Shift (MHz)'); xlabel('Distance (m)');
title('Detected Frequency Shift along the Fiber');

figure(2)
[Z_mesh, F_mesh] = meshgrid(z(1:Nz-M), lags_freq / 1e6);
surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none')
view(35, 45); colormap('jet'); colorbar;
xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)'); zlabel('Correlation');
title('3D Cross-Correlation Map');
rotate3d on;
% Distance
xlim([1400 1600]); 

% Frequency shift
ylim([-500 500]); 

% Correlation
zlim([-0.5 1]);

fprintf('--- Simulation successfully completed! ---\n');


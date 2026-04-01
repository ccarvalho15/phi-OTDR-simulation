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
M = round(d/dz);            % Number of points contained within one pulse width
alpha = 0;                  % Fiber attenuation (neglected for this simulation)

%% ------------ RANDOM REFRACTIVE INDEX PROFILE -----------
% Simulating the inhomogeneous nature of the fiber core using Gaussian noise
sigma_n = 2e-6;                % Standard deviation of refractive index fluctuations
delta_n = sigma_n*randn(1,Nz);
n = n_ave + delta_n;           % Final refractive index profile along the fiber

%% ------------ FRESNEL REFLECTION COEFFICIENT ------------
% Calculating local reflection coefficients (Rayleigh backscattering centers)
r = zeros(1,Nz);
for i = 1:Nz-1
    r(i) = (n(i)-n(i+1))/(n(i)+n(i+1));
end

%% ------------ PROPAGATION PHASE --------------------------
% Phase accumulation due to light propagation through the fiber
beta = 2*pi*n*nu0/c;        % Propagation constant
phi = cumsum(beta)*dz;      % Integrated phase along the fiber length

%% ------------ φOTDR TRACE GENERATION ---------------------
% Generating the complex electric field for the backscattered signal
E = zeros(1,Nz);

for k = 1:Nz-M
    sum_field = 0;
    for m = 0:M-1
        idx = k+m;
        % Vector summation of backscattered fields within the pulse volume
        sum_field = sum_field + r(idx)*exp(1j*2*phi(idx));
    end
    E(k) = exp(-alpha*z(k)/2)*sum_field;
end

%% ------------ STATISTICS CHECK ---------------------------
% Validating that the backscattered signal follows Rayleigh statistics
figure (1)
subplot(3,1,1)
amp = abs(E);
h1 = histogram(amp/mean(amp), 100,'Normalization','pdf');
hold on
x_amp = linspace(0, max(h1.BinEdges), 100);
rayleigh = (pi * x_amp/2) .* exp(-pi*x_amp.^2/4);
plot(x_amp, rayleigh, 'r', 'LineWidth', 2)
title('Amplitude Distribution (Rayleigh)')

subplot(3,1,2)
h2 = histogram(abs(E).^2/mean(abs(E).^2),100,'Normalization','pdf');
hold on
x_intensity = linspace(0, max(h2.BinEdges), 100);
y_exp = exp(-x_intensity);
plot(x_intensity, y_exp, 'r', 'Linewidth', 2)
title('Intensity Distribution (Exponential)')

subplot(3,1,3)
h3 = histogram(angle(E),100,'Normalization','pdf');
hold on
line([min(angle(E)) max(angle(E))], [1/(2*pi) 1/(2*pi)], 'Color', 'r', 'LineWidth', 2)
title('Phase Distribution (Uniform)')


%% SIMULATION OF STATIC MEASUREMENTS
freq_range = 1000e6;        % Frequency scanning range (1000 MHz)
delta_f = 10e6;             % Frequency step (10 MHz)
% Absolute optical frequency vector centered at nu0
f = nu0 + (-freq_range/2 : delta_f : freq_range/2);  
Nf = length(f);

% Pre-allocating matrices for reference and signal fields (Frequency x Distance)
E_ref = zeros(Nf, Nz-M);
E_sig = zeros(Nf, Nz-M);

% gamma = thermo-optic coefficient of the fiber
% eta = thermal expansian coefficient of the fiber
% epsilon = apllied strain
% delta_T = temperature variation

% Sensitivity coefficients for silica fiber at room temperature (300 K)
gamma = 9.1e-6;             % Thermo-optic coefficient
eta = 0.5e-6;               % Thermal expansion coefficient

% Temperature and strain profiles
delta_T = zeros(1, Nz);
epsilon = zeros(1, Nz);     

% Applying specific temperature variations at distinct locations
delta_T(round(2800/dz):round(2805/dz)) = -0.3; % 0.3 K decrease
delta_T(round(2950/dz):round(2960/dz)) = 0.2;  % 0.2 K increase

% Calculating refractive index changes due to temperature
delta_n_temp = (gamma + n_ave * eta) * delta_T;
n_temp = delta_n_temp + n;

% Calculating refractive index changes due to strain
delta_n_strain = n_ave * (1 - 0.1 * n_ave.^2) * epsilon;
n_strain = delta_n_strain + n;

%% ------------ FREQUENCY SCANNING LOOP ------------
% Simulating the backscattered field for each frequency in the scan
for f_idx = 1:Nf
    current_nu = f(f_idx);
    
    % Reference propagation phase (no temperature change)
    beta_ref = 2 * pi * n * current_nu /c;
    phi_ref = cumsum(beta_ref)*dz;      % Integrated phase

    % Signal propagation phase (including temperature effects)
    beta_sig = 2 * pi * n_temp * current_nu / c;
    phi_sig = cumsum(beta_sig)*dz;

    for k = 1:Nz-M
        idx = k:(k+M-1);
        % Complex field summation for both reference and measurement traces
        E_ref(f_idx, k) = sum(r(idx) .* exp(1j * 2 * phi_ref(idx)));
        E_sig(f_idx, k) = sum(r(idx) .* exp(1j * 2 * phi_sig(idx)));
    end
end

%% ------------ CROSS-CORRELATION ------------
% Calculating frequency shift via local cross-correlation of Rayleigh spectra
freq_shift = zeros(1, Nz-M); 
lags_freq = (-(Nf-1):(Nf-1)) * delta_f; % Frequency lag axis
corr_map = zeros(length(lags_freq), Nz-M); % Matrix for 3D visualization

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
figure(2)
[Z_mesh, F_mesh] = meshgrid(z(1:Nz-M), lags_freq / 1e6);
surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none')
view(35, 45); 
colormap('jet'); colorbar;
xlabel('Distance (m)'); ylabel('Frequency Shift (MHz)'); zlabel('Cross-correlation value (a.u)')
title('3D Cross-Correlation Map')
rotate3d on;

% Distance
xlim([2700 3000]); 

% Frequency shift
ylim([-500 500]); 

% Correlation
zlim([-0.5 1]);

fprintf('--- Simulation successfully completed! ---\n');


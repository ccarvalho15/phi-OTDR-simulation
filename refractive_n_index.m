clear; clc; close all

%% ----- 1. SYSTEM CONFIGURATION & WAVEGUIDE PROPERTIES -----
%%%%%%%%%%%%%%%%%%%%%%%%%
% - The fiber is treated as a series of inhomogeneities with random 
% refractive indices
%%%%%%%%%%%%%%%%%%%%%%%%%

% PHYSICAL AND OPTICAL CONSTANTS
c = 3e8;                    % Speed of light in vacuum (m/s)
lambda0 = 1550e-9;          % Operating wavelength (m)
nu0 = c/lambda0;            % Central optical frequency (Hz)

% FIBER PROPERTIES (modeled as a 1D waveguide)
n_ave = 1.456;              % Average refractive index of the silica fiber
L = 1000;                   % Total fiber length (m)
alpha = 0;                  % Fiber attenuation (neglected to explore intrinsic backscatter properties)

% SPATIAL SAMPLING AND RESOLUTION
dz = 0.05;                  % Spatial sampling interval (m)
z = 0:dz:L-dz;              % Distance vector along the fiber 1D model
Nz = length(z);             % Total number os spatial sampling points

% PULSE CHARACTERISTICS AND PHASE INTEGRATION
pulse_width = 10e-9;                % Temporal pulse width (10 ns)
d = c * pulse_width / (2 * n_ave);  % Spatial resolution (~1 m)
M = round(d / dz);                  % Number of scattering segments (inhomogeneities) within one pulse

fprintf('\n--- Fiber Simulation Initialization ---\n');
fprintf(['Fiber Length (L): %d m | ' ...
    'Spatial Resolution (d): %.2f m | ' ...
    'Total Spatial Points (Nz): %d | ' ...
    'Number of Scattering Segments (M): %d\n'], L, d, Nz, M);

%% ----- 2. GENERATION OF STOCHASTIC RAYLEIGH SCATTERING CENTERS -----
%%%%%%%%%%%%%%%%%%%%%%%%%
% - Modeling the fiber as a 1D waveguide with random inhomogeneities.
% - According to the paper, Rayleigh scattering is simulated by small 
% fluctuations in the refractive index along the fiber core.
%%%%%%%%%%%%%%%%%%%%%%%%%

sigma_n = 2e-6;                % Standard deviation of index fluctuations 
delta_n = sigma_n * randn(1, Nz); % Gaussian distributed random index variations
n = n_ave + delta_n;           % Resulting refractive index profile, n(z)

% FRESNEL REFLECTION COEFFICIENT
%%%%%%%%%%%%%%%%%%%%%%%%%
% - Calculation of the local reflection coefficient (r) at each interface.
% - The model treats each 'dz' step as a discrete boundary between media 
% with slightly different refractive indices.
%%%%%%%%%%%%%%%%%%%%%%%%%

r = zeros(1, Nz);
for i = 1:Nz-1
    % Fresnel formula for normal incidence between two adjacent cells
    % This represents the local backscattering amplitude at position z
    r(i) = (n(i) - n(i+1)) / (n(i) + n(i+1));
end

%% ----- 3. MODELING OF ENVIRONMENTAL SENSING EVENTS (STRAIN/TEMP) -----
%%%%%%%%%%%%%%%%%%%%%%%%%
% - Following the approach of Lu & Thomas, environmental perturbations 
% (e.g., strain or temperature) are modeled as localized modulations of the 
% refractive index profile.
%%%%%%%%%%%%%%%%%%%%%%%%%

delta_n_pert = zeros(1, Nz);

% CONFIGURATION OF SENSING EVENTS
num_events = 2;         % Number of discrete perturbation zones
pert_length = 5;        % Spatial width of each perturbation event (m)
spacing = 250;          % Spatial separation between events (m)

fprintf(['Number of events: %d | Width of each event: %.2f m | Spacing (Nz): ' ...
    '%d m\n'], num_events, pert_length, spacing);

% --- MAGNITUDE TRADE-OFF (delta_n vs. freq_range) ---
% PERTURBATION EQUATION: delta_nu / nu0 = - delta_n / n_ave
% The measured frequency shift (delta_nu) is directly proportional to index variation.
% RULE: A delta_n of 1e-4 causes shifts of ~13 GHz. To keep the signal within a 
% realistic observation window (e.g., 1 GHz) and avoid aliasing in the 
% cross-correlation algorithm, we use values ~10^-7 (tens of MHz).

min_mag = 1e-7; 
max_mag = 5e-7;

% Check if defined events fit within total fiber length L
total = (num_events * pert_length) + ((num_events - 1) * spacing);
if total > L
    error('Event configuration exceeds fiber length.')
end

% Randomized placement of the event sequence along the fiber
random_start = (L - total) * rand();
first_event = random_start;
last_event = random_start + total;

fprintf('\n--- Perturbation Events ---\n');
for i = 1:num_events
    % Mapping physical coordinates to vector indices
    start_idx = max(1, round(random_start/dz));
    end_idx = min(Nz, start_idx + round(pert_length/dz));

    random_end = random_start + pert_length;

    % Assign random magnitude and sign to simulate tensile or compressive stress
    random_mag = (min_mag + (max_mag - min_mag) * rand()) * sign(rand - 0.5);

    % Update the perturbation profile: delta_n(z)
    delta_n_pert(start_idx:end_idx) = random_mag;

    fprintf('Event %d :: Location: [%.2f; %.2f] m | Delta_n: %.2e\n', ...
        i, random_start, random_end, random_mag);

    % Increment position for the next event based on defined spacing
    random_start = random_start + pert_length + spacing;
end

% Final perturbed refractive index profile used for backscatter calculation
n_pert = n + delta_n_pert;

%% ----- 4. PROBE SIGNAL & FREQUENCY SWEEP PARAMETERS -----
%%%%%%%%%%%%%%%%%%%%%%%%%
% - In this stage, we simulate a frequency-swept probe signal to recover 
% the Rayleigh Backscatter (RB) spectra, as detailed in the static 
% measurement section of the paper.
%%%%%%%%%%%%%%%%%%%%%%%%%

freq_range = 1000e6;        % Total frequency scanning range (e.g., 1 GHz)
delta_f = 5e6;              % Frequency tuning step (5 MHz)

% Vector of absolute optical frequencies (nu) centered around the carrier frequency nu0
f = nu0 + (-freq_range/2 : delta_f : freq_range/2);  
Nf = length(f);

% Pre-allocation of matrices for the backscattered electric fields
% - Rows represent frequency components; 
% - Columns represent spatial positions (traces)
% The length is Nz-M+1 because the pulse integration window 'M' reduces the valid range.
E_ref = zeros(Nf, Nz-M+1);
E_sig = zeros(Nf, Nz-M+1);

% % The optical pulse is modeled as a rectangular window function of length M.
% This window determines the range of scattering centers that contribute 
% to the total interference at a given time delay.
window = ones(1, M); 

%% ----- 5. COHERENT BACKSCATTER INTEGRATION (PHASOR SUM) -----
%%%%%%%%%%%%%%%%%%%%%%%%%
% - This loop simulates the frequency-dependent backscattered field. 
% - For each frequency step, the total electric field is calculated by 
% integrating the contributions of all scattering centers within the pulse volume.
%%%%%%%%%%%%%%%%%%%%%%%%%

for f_idx = 1:Nf
    current_nu = f(f_idx);
    
    % --- Propagation Phase Calculation ---
    % beta represents the propagation constant. 
    % The phase phi is the spatial integral of beta along the fiber (cumsum).
    
    % Reference state: Cumulative phase with baseline refractive index n
    beta_ref = 2 * pi * n * current_nu /c;
    phi_ref = cumsum(beta_ref) * dz; 

    % Perturbed state: Cumulative phase with modified index n_pert
    %  This captures the phase shift induced by the sensing event (strain/temp).
    beta_sig = 2 * pi * n_pert * current_nu / c;
    phi_sig = cumsum(beta_sig) * dz;

    % --- Electric Field Construction via Convolution ---
    % WHY CONVOLUTION ('conv')? 
    % Physically, the detector captures the coherent phasor sum of all 
    % M reflectors within the pulse width at any given time delay. 
    % Mathematically, this is equivalent to a sliding window integration.
    
    % The term 'r .* exp(1j * 2 * phi)' represents the local backscattered 
    % light from each segment, where '2*phi' accounts for the round-trip path.
    
    % Using 'conv' with a rectangular 'window' replaces a nested spatial loop, 
    % significantly optimizing the simulation while maintaining exact 
    % physical consistency with the 1D waveguide model.
    
    E_ref(f_idx, :) = conv(r .* exp(1j * 2 * phi_ref), window, 'valid');
    E_sig(f_idx, :) = conv(r .* exp(1j * 2 * phi_sig), window, 'valid');
end

%% ----- 6. SPECTRAL SHIFT ESTIMATION VIA CROSS-CORRELATION -----
%%%%%%%%%%%%%%%%%%%%%%%%%
% - Following the methodology described in the paper, we perform a local 
% cross-correlation between the reference and perturbed Rayleigh backscatter 
% intensity spectra to quantify the environmental impact.
%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the frequency lag axis for correlation mapping
lags_freq = (-(Nf-1):(Nf-1)) * delta_f; 
corr_map = zeros(length(lags_freq), Nz-M+1); % Pre-allocate matrix for 3D correlation visualization
freq_shift = zeros(1, Nz-M+1); % Vector to store the estimated frequency shift per position

for k = 1:Nz-M+1
    % - Convert electric fields to intensity spectra (Power Spectral Density 
    % approximation)
    % - The paper analyzes the 'fading' pattern shifts in the frequency domain.
    ref = abs(E_ref(:,k)).^2;
    sig = abs(E_sig(:,k)).^2;
    
    % % --- Cross-Correlation Process ---
    % - We use zero-mean intensity spectra to eliminate DC bias and enhance 
    % the correlation peak detection.
    % - 'coeff' normalizes the sequences so that the auto-correlation at zero lag is 1.0.
    [cv, lags] = xcorr(sig - mean(sig), ref - mean(ref), 'coeff');
    corr_map(:, k) = cv;
    
    % --- Peak Tracking ---
    % The frequency shift (Delta_nu) corresponds to the lag that maximizes 
    % the correlation coefficient.
    [~, max_idx] = max(cv);
    freq_shift(k) = lags(max_idx) * delta_f; 
end

%% ----- 7. DATA VISUALIZATION & SENSOR PERFORMANCE ANALYSIS -----
%%%%%%%%%%%%%%%%%%%%%%%%%
% This section visualizes the mapping between the physical perturbation and 
% the recovered frequency shifts, simulating the output of a distributed 
% fiber sensing interrogation system.
%%%%%%%%%%%%%%%%%%%%%%%%%

% --- FIGURE 1: 1D Frequency Shift Profile ---
% This plot represents the demodulated sensing signal.
% The peaks should align with the 'start_idx' and 'end_idx' defined in Section 3.
figure(1)
plot(z(1:Nz-M+1), freq_shift / 1e6, 'LineWidth', 1.5)
grid on; ylabel('Frequency Shift (MHz)'); xlabel('Distance (m)');
title('Detected Frequency Shift along the Fiber');

%%
% --- FIGURE 2: 3D Correlation Surface ---
% Provides a global view of the cross-correlation peaks. 
% High correlation (close to 1.0) indicates high similarity between 
% the reference and signal Rayleigh spectra at the shifted frequency.
figure(2)
[Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none')
view(35, 45); colormap('jet'); colorbar;
xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)'); zlabel('Correlation');
title('3D Cross-Correlation Map');
rotate3d on;

%%
% --- FIGURE 3: 2D Contour Map with Peak Tracking ---
% A "top-down" view that combines the correlation energy with the 
% mathematical peak detection (white line).
% This visualization is excellent for assessing the Signal-to-Noise Ratio (SNR).
figure(3)
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
figure(4)
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

fprintf('\n--- Simulation successfully completed! ---\n');
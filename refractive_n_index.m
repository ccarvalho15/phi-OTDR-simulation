clear; clc; close all

%% ----- 1. SYSTEM CONFIGURATION & WAVEGUIDE PROPERTIES -----
clc; 
%%%%%%%%%%%%%%%%%%%%%%%%%
%   The fiber is treated as a series of inhomogeneities with random refractive 
% indices
%%%%%%%%%%%%%%%%%%%%%%%%%

% --- 1.1 PHYSICAL AND OPTICAL CONSTANTS ---
c = 3e8;                    % Speed of light in vacuum (m/s)
lambda0 = 1550e-9;          % Operating wavelength (m)
nu0 = c/lambda0;            % Central optical frequency (Hz)
n_ave = 1.456;              % Average refractive index of the silica fiber


% --- 1.2 FIBER PROPERTIES AND OPTICAL ATTENUATION (LOSS) SETUP ---
L = 240;                    % Total fiber length used on the lab (m) 

%   Attenuation is typically provided in dB/km (logarithmic scale). Standard 
% single-mode fiber (SMF-28) at 1550nm has approx. 0.19–0.2 dB/km range.
% References: 
%   1. D. A. A. Mello and F. A. Barbosa, Digital Coherent Optical Systems: 
% Architecture and Algorithms. Cham, Switzerland: Springer, 2021. doi: 
% 10.1007/978-3-030-66541-8.
%   2. R. Hui, Introduction to Fiber‑Optic Communications. London, UK: 
% Academic Press/Elsevier, 2020. ISBN: 978‑0‑12‑805345‑4.
attenuation = 0.2;

%   To use attenuation in the exponential field equations, we must convert 
% dB/km to the linear attenuation coefficient alpha (m^-1).
% References: 
%   1. G. Keiser, Optical Fiber Communications, 4th ed. New York, NY, USA: 
% McGraw‑Hill, 2011. ISBN: 978‑0‑07‑338071‑1.
%   2. R. Hui, Introduction to Fiber‑Optic Communications. London, UK: 
% Academic Press/Elsevier, 2020. ISBN: 978‑0‑12‑805345‑4.
alpha = attenuation/(10 * log10(exp(1)) * 1000);


% --- 1.3 SPATIAL SAMPLING AND RESOLUTION
dz = 0.05;                  % Spatial sampling interval (m)
z = 0:dz:L-dz;              % Distance vector along the fiber 1D model
Nz = length(z);             % Total number os spatial sampling points


% 1.4 --- PULSE AND MODULATION CHARACTERISTICS ---
pulse_width = 10e-9;                % Temporal pulse width (10 ns)
d = c * pulse_width / (2 * n_ave);  % Spatial resolution (~1 m)
M = round(d / dz);                  % Number of scattering segments 
                                    % (inhomogeneities) within one pulse

%   Coherent optical systems can work with roll-off factors ranging from 0.01 
% to 0.1, implementing a pulse with near-rectangular spectrum.
% References: 
%   1. D. A. A. Mello and F. A. Barbosa, Digital Coherent Optical Systems: 
% Architecture and Algorithms. Cham, Switzerland: Springer, 2021. doi: 
% 10.1007/978-3-030-66541-8.
roll_off = 0.1;                                      % Roll-off factor for pulse
                                                     % shaping (Raised Cosine)
t_rect = linspace(-pulse_width/2, pulse_width/2, M); % Time vector for pulse window


% --- 1.5 NOISE PARAMETERS ---
%   Laser linewidth: a narrower linewidth (e.g., 1 kHz) would increase 
% coherence; a wider one (> 1 MHz) would significantly increase phase noise 
% and degrade the correlation peak.
% References: 
%   1. S. Bottacchi, Noise and Signal Interference in Optical Fiber 
% Transmission Systems: An Optimum Design Approach. Hoboken, NJ, USA: 
% Wiley, 2009. ISBN: 978‑0‑470‑77056‑3.
linewidth = 100e3; % 100 kHz
%   SNR_dB: 20 dB is a common operating point for lab-bench systems, providing a 
% balance between sufficient signal strength and realistic noise levels to 
% test signal processing algorithms.
SNR_dB = 20; % in dB


fprintf('--- Fiber Simulation Initialization ---\n');
fprintf([ ...
    'Fiber Length (L):                      %d m\n' ...
    'Spatial Resolution (d):                %.2f m\n' ...
    'Total Spatial Points (Nz):             %d\n' ...
    'Number of Scattering Segments (M):     %d\n' ...
    'Attenuation (α):                       %.2f dB/km\n' ...
    'Roll-off factor (β):                   %.2f\n' ...
    'Laser Linewidth (Δν):                  %.1f kHz\n' ...
    'System Signal-to-Noise Ratio:          %d dB'], ...
    L, d, Nz, M, attenuation, roll_off, linewidth/1e3, SNR_dB);

%% ----- 2. GENERATION OF STOCHASTIC RAYLEIGH SCATTERING CENTERS -----
%%%%%%%%%%%%%%%%%%%%%%%%%
%   Modeling the fiber as a 1D waveguide with random inhomogeneities. According 
% to the paper, Rayleigh scattering is simulated by small fluctuations in 
% the refractive index along the fiber core.
%%%%%%%%%%%%%%%%%%%%%%%%%

sigma_n = 2e-6;                     % Standard deviation of index fluctuations 
delta_n = sigma_n * randn(1, Nz);   % Gaussian distributed random index variations
n = n_ave + delta_n;                % Resulting refractive index profile, n(z)

% FRESNEL REFLECTION COEFFICIENT
%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculation of the local reflection coefficient (r) at each interface.
%   The model treats each 'dz' step as a discrete boundary between media 
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
%   Following the approach of Lu & Thomas, environmental perturbations 
% (e.g., strain or temperature) are modeled as localized modulations of the 
% refractive index profile.
%   Based on Lu & Thomas, temperature changes modulate the phase via dn/dT.
% For Silica: dn/dT approx. 1.1e-5 
%%%%%%%%%%%%%%%%%%%%%%%%%

delta_n_pert = zeros(1, Nz);

% CONFIGURATION OF SENSING EVENTS
num_events = 5;         % Number of discrete perturbation zones
pert_length = 1.5;        % Spatial width of each perturbation event (m)
spacing = 4;            % Spatial separation between events (m)

fprintf([ ...
    '\nNumber of events:                      %d\n' ...
    'Width of each event:                   %.2f m\n' ...
    'Spacing:                               %.2f m\n'], ...
    num_events, pert_length, spacing);

%%
% --- MAGNITUDE TRADE-OFF (delta_n vs. Correlation Bandwidth) ---
%   The speckle pattern generates a correlation peak with a finite bandwidth
% (FWHM).
%   For a pulse width of 10 ns (d ~ 1.03 m), the bandwidth is approximately:
% FWHM = c / (2 * n_ave * d) ~= 100 MHz.
%
% Option A (Realistic): delta_n = 1e-7 -> Shift ~= 19 MHz.
% The shift occurs *within* the correlation peak width, appearing as a slight 
% deformation in the 3D surface plot.
%
% Option B (Visual): delta_n > 7.5e-7 -> Shift > 100 MHz.
% The shift moves the peak *outside* the original bandwidth, making the 
% sensing event visually distinct in the 'surf' and 'contour' maps.

% --- INTERACTIVE MAGNITUDE SELECTION (CONSOLE) ---
fprintf('\nSelection of Perturbation Magnitude:         1) Realistic (Micro-events)        2) Visual (Distinct Peaks)\n');
user_choice = input('> Option ');

if isempty(user_choice) || user_choice == 1
    min_mag = 1e-7; 
    max_mag = 5e-7;
elseif user_choice == 2
    min_mag = 5e-7; 
    max_mag = 1e-6;
else
    min_mag = 1e-7; 
    max_mag = 5e-7;
    fprintf('Invalid input. Defaulting to Realistic Micro-events.\n');
end
%% 
% Defining the sensing zone starting point (e.g., L-30 meters)
sensing_zone = L - 30;

% Check if defined events fit within total fiber length L
total = (num_events * pert_length) + ((num_events - 1) * spacing);
if total > 30
    error(['Event configuration exceeds the allocated 30m sensing zone. Your ' ...
        'current setup needs %.2f m.'], total);
end

% Randomized placement of the event sequence along the fiber
random_start = sensing_zone + (30 - total) * rand();
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
%   In this stage, we simulate a frequency-swept probe signal to recover 
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
% The length is Nz-M+1 because the pulse integration window 'M' reduces the 
% valid range.
E_ref = zeros(Nf, Nz-M+1);
E_sig = zeros(Nf, Nz-M+1);

% To capture the ideal electric field
E_ref_id = zeros(Nf, Nz-M+1); 
E_sig_id = zeros(Nf, Nz-M+1); 

%% ----- 4.1 REALISTIC PULSE SHAPPING -----
%   Instead of a perfecdt retangular pulse, using window  = ones(1, M) where 
% the optical pulse is modeled as a rectangular window function of length M,
% we use a Raised Cosine window to model the Electro-Optic Modulator (EOM) 
% rise and fall times
%
% The time-domain Raised Cosine pulse models the EOM output envelope. It is 
% the inverse Fourier transform of the raised cosine frequency response:
% h(t) = sinc(t/T) * cos(pi * roll_off * t / pulse_width) / 
%                                          (1 - (2*roll_off*t/pulse_width)^2)
% 

t_norm = t_rect / pulse_width; % Normalized time
sinc_term = sinc(t_norm); % sinc time
cos_term = cos(pi * roll_off * t_norm); % cosine modulation term
denom = 1 - (2 * roll_off * t_norm).^2; % denominator: goes to zero at t = +-T/(2*roll_off)

% Compute the raised cosine pulse h(t)
window = sinc_term .* cos_term ./ denom;

% Using L' Hospital's rule to handle the two singularities
% lim_{t -> t+- T/(2*beta} h(t) = (pi/4) * sinc(1/2*beta))
singularity_mask = abs(denom) < 1e-6;
window(singularity_mask) = (pi/4) * sinc(1 / (2 * roll_off));

% Normalize the amplitude to 1 to maintain consistency in backscatter intensity
window = window / max(abs(window));

%% ----- 5. COHERENT BACKSCATTER INTEGRATION & NOISE MODELING -----
%%%%%%%%%%%%%%%%%%%%%%%%%
%   This loop simulates the frequency-dependent backscattered field. For 
% each frequency step, the total electric field is calculated by integrating 
% the contributions of all scattering centers within the pulse volume.
%
% NOISE MODEL:
% - Laser Phase Noise: Modeled as a Wiener Process (aka Random Walk), the
% phase uncertainty accumulates with time (and thus distance).
% - Receiver Noise: Combines thermal and shot noise, and it's added as
% complex additive white Guassian noise to the field
%%%%%%%%%%%%%%%%%%%%%%%%%

for f_idx = 1:Nf
    current_nu = f(f_idx);
    

    % --- 5.1 Propagation Phase Calculation ---
    % beta: propagation constant
    % phi: spatial integral of beta along the fiber (cumsum).
    
    % Reference state: Cumulative phase with baseline refractive index n
    beta_ref = 2 * pi * n * current_nu /c;
    phi_ref = cumsum(beta_ref) * dz; 

    % Perturbed state: Cumulative phase with modified index n_pert
    % This captures the phase shift induced by the sensing event (strain/temp).
    beta_sig = 2 * pi * n_pert * current_nu / c;
    phi_sig = cumsum(beta_sig) * dz;


    % --- 5.2 Electric Field Construction via Convolution ---
    % WHY CONVOLUTION ('conv')? 
    % Physically, the detector captures the coherent phasor sum of all 
    % M reflectors within the pulse width at any given time delay. 
    % Mathematically, this is equivalent to a sliding window integration.
    
    % The term 'r .* exp(1j * 2 * phi)' represents the local backscattered 
    % light from each segment, where '2*phi' accounts for the round-trip path.
    
    % Using 'conv' with a rectangular 'window' replaces a nested spatial loop, 
    % significantly optimizing the simulation while maintaining exact 
    % physical consistency with the 1D waveguide model.

    E_ref_conv = conv(r .* exp(1j * 2 * phi_ref), window, 'valid');
    E_sig_conv = conv(r .* exp(1j * 2 * phi_sig), window, 'valid');

    % Adjust distance vector for 'valid' conv
    z_valid = z(1:length(E_ref_conv)); 


    % --- 5.3 Laser Phase Noise (Transmitter Impairment) ---
    % Phase noise variance increases linearly with the round-trip delay (tau) 
    tau = 2 * n_ave * z_valid / c; % Time delay
    std_dvt_phase = sqrt(2 * pi * linewidth * tau); % sqrt(2 * pi * linewidth * tau)
    
    % Generate independent phase noise for this specific frequency step
    laser_phase_noise_ref = std_dvt_phase .* randn(size(E_ref_conv));
    laser_phase_noise_sig = std_dvt_phase .* randn(size(E_sig_conv));
        
    % Apply the stochastic phase jitter to the complex fields
    E_ref_conv = E_ref_conv .* exp(1j * laser_phase_noise_ref);
    E_sig_conv = E_sig_conv .* exp(1j * laser_phase_noise_sig);

    
    % --- 5.4 Fiber Loss (Beer-Lambert Law) ---
    % As the light travels to distance 'z' and back to the detector (round-trip),
    % the electric field amplitude decays exponentially.
    % We use 'exp(-alpha * z)' because the signal accumulates loss over the 
    % total path (2*z). Since alpha is defined for power, the field decay 
    % over distance '2z' is exp(-(alpha/2) * 2z) = exp(-alpha * z).
    loss_factor = exp(-alpha * z_valid);

    E_ref_raw = loss_factor .* E_ref_conv;
    E_sig_raw = loss_factor .* E_sig_conv;


    % --- 5.5 Additive White Gaussian Noise (Receiver Impairment) ---
    % Simulates the electronic noise floor (Thermal) and photon counting 
    % noise (Shot). The 'measured' flag ensures the noise power is scaled 
    % relative to the signal power.
    E_ref(f_idx, :) = awgn(E_ref_raw, SNR_dB, 'measured');
    E_sig(f_idx, :) = awgn(E_sig_raw, SNR_dB, 'measured');

    % Store ideal fields for performance benchmarking
    E_ref_id(f_idx, :) = E_ref_conv;
    E_sig_id(f_idx, :) = E_sig_conv;

end

%% ----- 6. SPECTRAL SHIFT ESTIMATION VIA CROSS-CORRELATION -----
%%%%%%%%%%%%%%%%%%%%%%%%%
%   Following the methodology described in the paper, we perform a local 
% cross-correlation between the reference and perturbed Rayleigh backscatter 
% intensity spectra to quantify the environmental impact.
%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the frequency lag axis for correlation mapping
lags_freq = (-(Nf-1):(Nf-1)) * delta_f; 
corr_map = zeros(length(lags_freq), Nz-M+1); % Pre-allocate matrix for 3D correlation visualization
freq_shift = zeros(1, Nz-M+1);               % Vector to store the estimated frequency shift per position

for k = 1:Nz-M+1
    % - Convert electric fields to intensity spectra (Power Spectral Density 
    % approximation)
    % - The paper analyzes the 'fading' pattern shifts in the frequency domain.
    ref = abs(E_ref(:,k)).^2;
    sig = abs(E_sig(:,k)).^2;
    
    % % --- Cross-Correlation Process ---
    %   We use zero-mean intensity spectra to eliminate DC bias and enhance 
    % the correlation peak detection.
    %   'coeff' normalizes the sequences so that the auto-correlation at zero lag is 1.0.
    [cv, lags] = xcorr(sig - mean(sig), ref - mean(ref), 'coeff');
    corr_map(:, k) = cv;
    
    % --- Peak Tracking ---
    %   The frequency shift (Delta_nu) corresponds to the lag that maximizes 
    % the correlation coefficient.
    [~, max_idx] = max(cv);
    freq_shift(k) = lags(max_idx) * delta_f; 
end

%% ----- 7. DATA VISUALIZATION & SENSOR PERFORMANCE ANALYSIS -----
%%%%%%%%%%%%%%%%%%%%%%%%%
%   This section visualizes the mapping between the physical perturbation and 
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
xlim([max(210, first_event - 10) min(L, last_event + 10)]);  % Fiber length
ylim([-250 250]); % Frequency shift window
zlim([-0.5 1]); % Correlation magnitude scale

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

%% ----- 8. ATTENUATION IMPACT ANALYSIS -----

% Power configuration
P_input_mW = 10; % Source power in milliwatts (10 mW)
% To convert mW to dBm: P [dBm] = 10 * log10 (P [mW] / 1 [mW])
P_input_dBm = 10 * log10(P_input_mW / 1); % Logarithmic power reference (10 dBm)

% Normalizing electric field
% Calculates the scaling factor to map dimensionless simulated fields to 
% physical mW. Using mean of the first point to stabilize against coherent 
% fading (speckle).
scale_factor = P_input_mW / mean(abs(E_ref_id(1,:))).^2;

% Select the first frequency from the sweep for visualization
P_ref_ideal = abs(E_ref_id(1, :)).^2 * scale_factor;
P_ref_loss = abs(E_ref(1, :)).^2 * scale_factor;

% Comparing the backscattered intensity with and without fiber loss
figure (5);
z_axis = z(1:length(P_ref_loss));

% Subplot 1: Linear Scale
subplot(2,1,1);
plot(z_axis, P_ref_ideal, 'b', 'DisplayName', 'Ideal (No Loss)');
hold on;
plot(z_axis, P_ref_loss, 'r', 'DisplayName', ['Fiber Loss (', num2str(attenuation), ' dB/km) + Noise']);
hold off;
title('Backscattered Intensity (Linear Scale)'); 
xlabel('Distance (m)'); ylabel('Power (mW)'); legend('Location', 'northeast'); grid on;

% Subplot 2: Logarithmic Scale (OTDR Trace)
subplot(2,1,2);
plot(z_axis, 10*log10(P_ref_ideal + eps), 'b', 'DisplayName', 'Ideal (No Loss)');
hold on;
plot(z_axis, 10*log10(P_ref_loss + eps), 'r', 'DisplayName', ['Fiber Loss (', num2str(attenuation), ' dB/km) + Noise']);
hold off;
title('Backscattered Intensity (Logarithmic Scale)'); 
ylabel('Power (dBm)'); xlabel('Distance (m)'); grid on;  legend('Location', 'northeast');

fprintf('\n--- Simulation successfully completed! ---\n');


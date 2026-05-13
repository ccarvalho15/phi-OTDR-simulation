clear; clc; close all

%% ========================================================================
% 1. SYSTEM CONFIGURATION & WAVEGUIDE PROPERTIES
% ========================================================================
clc; 
%   The fiber is treated as a series of inhomogeneities with random refractive 
% indices


% --- 1.1 PHYSICAL AND OPTICAL CONSTANTS ---
c = 3e8;                    % Speed of light in vacuum (m/s)
lambda0 = 1550e-9;          % Operating wavelength (m)
nu0 = c/lambda0;            % Central optical frequency (Hz)


% --- 1.2 FIBER MEDIUM PROPERTIES ---
L = 240;           % Total fiber length used on the lab (m) 
n_ave = 1.456;     % Average refractive index of the silica fiber
attenuation = 0.2; % Attenuation is typically provided in dB/km (logarithmic 
                   % scale). Standard single-mode fiber (SMF-28) at 1550nm 
                   % has approx. 0.19–0.2 dB/km range.
% References: [1] Digital Coherent Optical Systems: Architecture and Algorithms. 2021
%             [2] Introduction to Fiber‑Optic Communications. 2020.


% --- 1.3 LOSS & POWER CALIBRATION ---
%   To use attenuation in the exponential field equations, we must convert 
% dB/km to the linear attenuation coefficient alpha (m^-1).
alpha = attenuation/(10 * log10(exp(1)) * 1000); 
% References: [1] Optical Fiber Communications. 2011.
%             [2] Introduction to Fiber‑Optic Communications. 2020.

%   To convert mW to dBm: P [dBm] = 10 * log10 (P [mW] / 1 [mW])
P_input_mW = 10;                          % Source power in milliwatts (10 mW)
P_input_dBm = 10 * log10(P_input_mW / 1); % Logarithmic power reference (10 dBm)


% --- 1.4 DISCRETE SAMPLING ---
dz = 0.05;       % Spatial sampling interval (m)
z = 0:dz:L-dz;   % Distance vector along the fiber 1D model
Nz = length(z);  % Total number os spatial sampling points


% --- 1.5 PULSE-LIMITED RESOLUTION ---
pulse_width = 10e-9;                % Temporal pulse width (10 ns)
d = c * pulse_width / (2 * n_ave);  % Spatial resolution (~1 m)
M = round(d / dz);                  % Number of scattering segments 
                                    % (inhomogeneities) within one pulse


% --- 1.6 PHASE & AMPLITUDE NOISE ---
linewidth = 100e3;     % Laser linewidth (100 kHz)
SNR_dB = 20;           % System Signal-to-Noise Ratio (dB)
sigma_n = 2e-6;        % Standard deviation of index fluctuations 

%   Laser linewidth: a narrower linewidth (e.g., 1 kHz) would increase 
% coherence; a wider one (> 1 MHz) would significantly increase phase noise 
% and degrade the correlation peak.
% References: [1] Noise and Signal Interference in Optical Fiber Transmission Systems: 
%                 An Optimum Design Approach. 2009.
%   SNR_dB: 20 dB is a common operating point for lab-bench systems, providing a 
% balance between sufficient signal strength and realistic noise levels to 
% test signal processing algorithms.


% --- 1.7 FREQUENCY SWEEP (SPECTRAL SCAN) ---
freq_range = 1000e6;        % Total frequency scanning range (e.g., 1 GHz)
delta_f = 5e6;              % Frequency tuning step (5 MHz)
% Vector of absolute optical frequencies (nu) centered around the carrier frequency nu0
f = nu0 + (-freq_range/2 : delta_f : freq_range/2);  
Nf = length(f);


% --- 1.8 SPATIAL DISTRIBUTION ---
num_events  = 5;        % Número de zonas de perturbação
pert_length = 1.5;      % Spatial width of each perturbation event (m)
spacing = 4;            % Spatial separation between events (m)
sensing_zone = L - 30;  % Defining the sensing zone starting point (e.g., 
                        % L-30 meters)

fprintf(' ==== PHI-OTDR SIMULATION CONFIGURATION ==== \n');
fprintf([ ...
    'Fiber Length (L):                      %d m                           Input Power:                           %d mW (%.1f dBm)\n' ...
    'Spatial Resolution (d):                %.2f m                          Sweep Range:                           %d MHz\n' ...
    'Total Spatial Points (Nz)              %d                            Frequency step (Δf):                   %d MHz\n' ...
    'Number of Scattering Segments (M):     %d                              Number of frequencies (Nf):            %d\n' ...
    'Sampling Interval (dz):                %.3f m                         Laser Linewidth (Δν):                  %.1f kHz\n' ...                   
    'Attenuation:                           %.2f dB/km                      System SNR:                            %d dB\n' ...
    'Attenuation Coefficient (α):           %.4e m^-1                 Index standard deviation (σ_n):        %.1e\n' ...
    'Average Refractive Index (n_ave):      %.3f                           Number of Events:                      %d\n' ...
    'Operating Wavelength (lambda_0):       %d nm                         Width of each event:                   %.2f m\n'...
    'Central frequency (nu_0)               %.4e Hz                   Separation between events:             %.2f m\n'...
    'Pulse Width:                           %d ns                           Sensing Zone Start:                    %d m\n'], ...
    L, P_input_mW, P_input_dBm, d, freq_range/1e6, Nz, delta_f/1e6, M, Nf, ...
    dz, linewidth/1e3, attenuation, SNR_dB, alpha, sigma_n, n_ave, num_events, ...
    lambda0*1e9, pert_length, nu0, spacing, pulse_width*1e9, sensing_zone);

        

%% ========================================================================
% 2. GENERATION OF STOCHASTIC RAYLEIGH SCATTERING CENTERS 
% ========================================================================
%   Modeling the fiber as a 1D waveguide with random inhomogeneities. According 
% to the paper, Rayleigh scattering is simulated by small fluctuations in 
% the refractive index along the fiber core.


delta_n = sigma_n * randn(1, Nz);   % Gaussian distributed random index variations
n = n_ave + delta_n;                % Resulting refractive index profile, n(z)

%% ========================================================================
% 2.1 FRESNEL REFLECTION COEFFICIENT
% ========================================================================
%   Calculation of the local reflection coefficient (r) at each interface.
%   The model treats each 'dz' step as a discrete boundary between media 
% with slightly different refractive indices.


r = zeros(1, Nz);
for i = 1:Nz-1
    % Fresnel formula for normal incidence between two adjacent cells
    % This represents the local backscattering amplitude at position z
    r(i) = (n(i) - n(i+1)) / (n(i) + n(i+1));
end

%% ========================================================================
% 3. MODELING OF ENVIRONMENTAL SENSING EVENTS (STRAIN/TEMP)
% ========================================================================
%   Following the approach of Lu & Thomas, environmental perturbations 
% (e.g., strain or temperature) are modeled as localized modulations of the 
% refractive index profile.
%   Based on Lu & Thomas, temperature changes modulate the phase via dn/dT.
% For Silica: dn/dT approx. 1.1e-5 

delta_n_pert = zeros(1, Nz);

magnitudes = [1e-7; -3e-7; 2e-7; 5e-7; -4e-7];

first_event = sensing_zone; 
last_event  = sensing_zone + (num_events-1)*(pert_length + spacing) + pert_length;

fprintf('\n--- Perturbation Events ---\n');
fprintf('\n%-8s | %-18s | %-12s | %-12s\n', 'Event', 'Location (m)', 'Delta_n', 'Shift (MHz)');
fprintf('%s\n', repmat('-', 1, 60));

for i = 1:num_events
    start_pos = 2 + sensing_zone + (i-1) * (pert_length + spacing);
    end_pos = start_pos + pert_length;
    mag = magnitudes(i);

    shift_MHz = (nu0 * mag / n_ave) / 1e6;

    % Mapping physical coordinates to vector indices
    start_idx = max(1, round(start_pos / dz));
    end_idx = min(Nz, round(end_pos / dz));

    % Update the perturbation profile: delta_n(z)
    delta_n_pert(start_idx:end_idx) = mag;

    loc_str = sprintf('[%.2f; %.2f]', start_pos, end_pos);
    fprintf('Event %-2d | %-18s | %-12.2e | %-12.2f\n', ...
        i, loc_str, mag, shift_MHz);

end

% Final perturbed refractive index profile used for backscatter calculation
n_pert = n + delta_n_pert;
fprintf('%s\n', repmat('-', 1, 60));

%% ========================================================================
% 4. PROBE SIGNAL & FREQUENCY SWEEP PARAMETERS 
% ========================================================================
%   In this stage, we simulate a frequency-swept probe signal to recover 
% the Rayleigh Backscatter (RB) spectra, as detailed in the static 
% measurement section of the paper.


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

% To analyze SNR between ideal signal and degraded signal
E_ref_raw_all = zeros(Nf, Nz-M+1);
E_sig_raw_all = zeros(Nf, Nz-M+1);

%% ========================================================================
% 4.1 REALISTIC PULSE SHAPING 
% ========================================================================

fprintf('\nSelection of Pulse Shapping:         1) RC Filter Model        2) Super-Gaussian Model\n');
shape_choice = input('> Option ');
if shape_choice == 1
    % ----- 4.1.a REALISTIC PULSE SHAPING (RC FILTER MODEL) -----
    % This section models the Electro-Optic Modulator (EOM) response as a 
    % first-order RC low-pass filter. This simulates the hardware constraints 
    % where the electronic driver cannot switch instantaneously, resulting 
    % in finite rise and fall times.
    
    % Filter configuration
    fc = 100e6;                     % Filter cutoff frequency (100 MHz)
    RC = 1 / (2 * pi * fc);         % RC time constant (seconds)
    tau_pulse = pulse_width;        % Duration of the ideal rectangular pulse (τ)
    A = 1;                          % Normalized peak amplitude
    
    % Time vector: Covers both the charging (rise) and discharging (fall) phases.
    % Although the physical pulse duration is tau_pulse, the RC tail extends 
    % beyond it. We generate 2*M points to visualize the full waveform.
    t_full = linspace(0, 2 * pulse_width, 2 * M);
    
    % Reference: Ideal rectangular input for comparison
    rect_pulse = ones(1, M);
    
    % Initialize the full waveform vector
    window_full = zeros(1, 2 * M);
    
    % Apply the step-response equations for an RC circuit
    for k = 1:2*M
        tk = t_full(k);
        if tk <= 0
            window_full(k) = 0;
        elseif tk < tau_pulse
            % Rising Edge (Charging Phase): 0 < t < τ
            % Equation: V(t) = A * (1 - exp(-t/RC))
            window_full(k) = A * (1 - exp(-tk / RC));
        else
            % Falling Edge (Discharging Phase): t > τ
            % The decay starts from the amplitude reached at the end of the pulse.
            % Equation: V(t) = V(tau) * exp(-(t - τ)/RC)
            window_full(k) = A * (1 - exp(-tau_pulse / RC)) * exp(-(tk - tau_pulse) / RC);
        end
    end
    
    % The 'window' used for the backscatter convolution represents the pulse 
    % within the sampling window M. 
    window = window_full(1:M);
    
    % Normalize the window to ensure peak power is consistent across models
    window = window / max(window);  

    % ---- VISUAL VERIFICATION (RC PULSE MODEL) ----
    figure(1)
    set(gcf, 'Name', 'Pulse Shape')
    
    % Time vector for display purposes only (covering twice the pulse width)
    t_plot = linspace(0, 2 * pulse_width * 1e9, 2 * M);  % Time in ns
    t_window_ns = linspace(0, pulse_width * 1e9, M);    % Time in ns for the active window
    
    
    % Shows the full charging/discharging cycle of the RC filter
    plot(t_plot, window_full, 'Color', [0.1 0.2 0.8], 'LineWidth', 2, ...
     'DisplayName', 'RC Pulse Shapping');
    hold on;
    plot(t_window_ns, ones(1,M), 'r--', 'LineWidth', 1.5, 'DisplayName', ...
        'Ideal Rectangular Pulse');
    hold off;

    xlabel('Time (ns)'); ylabel('Amplitude');
    title('RC Pulse Shapping (fc = 100 MHz)');
    legend('Location', 'northeast');
    ylim([0 1.2]); grid on;
    
elseif shape_choice == 2
    % ----- 4.1.b REALISTIC PULSE SHAPING (SUPER-GAUSSIAN MODEL) -----
    % The Super-Gaussian function is used to model pulses that have a 
    % smoother transition than a hard rectangular pulse but can maintain 
    % a "flat-top" characteristic depending on the order (N).
    
    % Pulse Order (N):
    % N = 1 results in a standard Gaussian shape.
    % N > 3 approaches a rectangular "flat-top" pulse with realistic edges.
    order_N = 3; 
    
    % Time Vector: Centered at zero to satisfy the Super-Gaussian symmetry.
    % It spans the full pulse duration defined by 'pulse_width'.
    t = linspace(-pulse_width, pulse_width, M);
   
    % Generate the Super-Gaussian window using the sgauss function.
    % Parameters: (time, FWHM width, Energy, Chirp, Order)
    % Here, Tfwhm is set to pulse_width to define the effective spatial cell.
    window = sgauss(t, pulse_width, 1, 0, order_N);
    
    % Normalization: Ensures the peak power is 1.0 to maintain 
    % consistency in the backscattered signal intensity across different models.
    window = window / max(window);
    
    % Reference: Create a standard rectangular pulse for visual comparison.
    rect_pulse = ones(1, M); 

    % ---- VISUAL VERIFICATION (SUPER-GAUSSIAN SHAPE) ----
    figure(1)
    set(gcf, 'Name', 'Pulse Shape')
    
    % Time vector in nanoseconds for plotting
    t_ns = t * 1e9;
    
    plot(t_ns, window, ...
        'Color', [0.2 0.2 0.8], 'LineWidth', 2, ...
        'DisplayName', ['Super-Gaussian (N=' num2str(order_N) ')']);
    hold on;
    plot(t_ns, rect_pulse, 'r--', 'LineWidth', 1.5, ...
         'DisplayName', 'Ideal Rectangular Pulse');
    hold off;
    title('Super Gaussian Model Pulse Shapping');
    xlabel('Relative Time (ns)'); ylabel('Normalized Amplitude');
    legend('Location', 'northeast');
    ylim([0 1.2]); 
    grid on;
end

%% ========================================================================
% 5. COHERENT BACKSCATTER INTEGRATION & NOISE MODELING
% ========================================================================
%   This loop simulates the frequency-dependent backscattered field. For 
% each frequency step, the total electric field is calculated by integrating 
% the contributions of all scattering centers within the pulse volume.
%
% NOISE MODEL:
% - Laser Phase Noise: Modeled as a Wiener Process (aka Random Walk), the
% phase uncertainty accumulates with time (and thus distance).
% - Receiver Noise: Combines thermal and shot noise, and it's added as
% complex additive white Guassian noise to the field

z_valid = z(1 : (Nz - M + 1));
t_laser = (0:Nz-M) * (2 * n_ave * dz / c);  % passo = round-trip time por célula

for f_idx = 1:Nf
    current_nu = f(f_idx);
    shift = current_nu - nu0;
    
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

    % --- 5.3 Laser Phase Noise (Transmitter Impairment) ---
    % Phase noise variance increases linearly with the round-trip delay (tau) 
    E_laser = lasercw(t_laser, 10, 0, linewidth, shift);

    % --- 5.4 Fiber Loss (Beer-Lambert Law) ---
    % As the light travels to distance 'z' and back to the detector (round-trip),
    % the electric field amplitude decays exponentially.
    % We use 'exp(-alpha * z)' because the signal accumulates loss over the 
    % total path (2*z). Since alpha is defined for power, the field decay 
    % over distance '2z' is exp(-(alpha/2) * 2z) = exp(-alpha * z).
    loss_factor = exp(-alpha * z_valid);

    E_ref_raw = E_ref_conv .* E_laser .* loss_factor;
    E_sig_raw = E_sig_conv .* E_laser .* loss_factor;

    % --- 5.5 Additive White Gaussian Noise (Receiver Impairment) ---
    % Simulates the electronic noise floor (Thermal) and photon counting 
    % noise (Shot). The 'measured' flag ensures the noise power is scaled 
    % relative to the signal power.
    E_ref(f_idx, :) = awgn(E_ref_raw, SNR_dB, 'measured');
    E_sig(f_idx, :) = awgn(E_sig_raw, SNR_dB, 'measured');

    % Store ideal fields for performance benchmarking
    E_ref_id(f_idx, :) = E_ref_conv;
    E_sig_id(f_idx, :) = E_sig_conv;

    % Store no ideal fields for performance benchmarking
    E_ref_raw_all(f_idx, :) = E_ref_raw;
    E_sig_raw_all(f_idx, :) = E_sig_raw;

end

%% ========================================================================
% 6. SPECTRAL SHIFT ESTIMATION VIA CROSS-CORRELATION 
% ========================================================================
%   Following the methodology described in the paper, we perform a local 
% cross-correlation between the reference and perturbed Rayleigh backscatter 
% intensity spectra to quantify the environmental impact.


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


%% ========================================================================
%  7. DATA VISUALIZATION & SENSOR PERFORMANCE ANALYSIS -----
% ========================================================================
%   This section visualizes the mapping between the physical perturbation and 
% the recovered frequency shifts, simulating the output of a distributed 
% fiber sensing interrogation system.


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

figure(3)
set(gcf,  'Name', '3D CC (reduced)');
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

P_ideal_mean = mean(abs(E_ref_id).^2, 1);
P_noisy_mean = mean(abs(E_ref).^2,    1);

% Normalizing electric field
% Calculates the scaling factor to map dimensionless simulated fields to 
% physical mW. Using mean of the first point to stabilize against coherent 
% fading (speckle).
scale_factor = P_input_mW / max(P_ideal_mean);

% Select the first frequency from the sweep for visualization
P_ref_ideal = P_ideal_mean * scale_factor;
P_ref_loss  = P_noisy_mean * scale_factor;

% Comparing the backscattered intensity with and without fiber loss
figure(6)
set(gcf, 'Name', 'Attenuation Impact');

% Subplot 2: Logarithmic Scale (OTDR Trace)
plot(z_valid, 10*log10(P_ref_ideal + eps), 'b', 'DisplayName', 'Ideal (No Loss)');
hold on;
plot(z_valid, 10*log10(P_ref_loss + eps), 'r', 'DisplayName', ['Fiber Loss (', num2str(attenuation), ' dB/km) + Noise']);
hold off;
title('Backscattered Intensity (Logarithmic Scale)'); 
ylabel('Power (dBm)'); xlabel('Distance (m)'); grid on;  legend('Location', 'southeast');


%%

fprintf('\n--- Simulation successfully completed! ---\n');


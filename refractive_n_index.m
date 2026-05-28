clear; clc; close all
addpath('functions\')

%% =======================================================================
% 1. SYSTEM CONFIGURATION & WAVEGUIDE PROPERTIES
% ========================================================================
%
%   The fiber is modelled as a 1D waveguide with random refractive-index
%   inhomogeneities. Each spatial step 'dz' represents a discrete boundary
%   between adjacent media, whose backscatter is computed via the Fresnel
%   formula. Environmental events (strain / temperature) are superimposed
%   as localised perturbations on the baseline index profile.
%
%   Key references:
%     [1] X. Lu & P. J. Thomas, "Numerical modeling of ΦOTDR sensing using
%         a refractive index perturbation approach", 2020.
%     [2] Digital Coherent Optical Systems: Architecture and Algorithms, 2021.
%     [3] Optical Fiber Communications, 2011.
%     [4] Introduction to Fiber-Optic Communications, 2020.
%     [5] Noise and Signal Interference in Optical Fiber Transmission
%         Systems: An Optimum Design Approach, 2009.


% ------------------------------------------------------------------------
% 1.1  FUNDAMENTAL PHYSICAL CONSTANTS
% ------------------------------------------------------------------------
c = 3e8;                    % Speed of light in vacuum (m/s)
lambda0 = 1550e-9;          % Operating wavelength (m)
nu0 = c/lambda0;            % Central optical frequency (Hz)


% -------------------------------------------------------------------------
% 1.2  FIBER MEDIUM PROPERTIES
% -------------------------------------------------------------------------
L = 240;           % Total fiber length used in the lab (m)
n_ave = 1.456;     % Average refractive index of the silica core
gamma = 9.1e-6;    % Thermo-optic coefficient  dn/dT (K^-1)
eta = 0.5e-6;      % Thermal expansion coefficient (K^-1)

%    Attenuation is typically provided in dB/km (logarithmic scale). Standard 
%    single-mode fiber (SMF-28) at 1550nm has approx. 0.19–0.2 dB/km range. 
%    [Refs 3,4]
attenuation = 0.2; % Fiber attenuation (dB/km)


% -------------------------------------------------------------------------
% 1.3  LOSS CONVERSION & INPUT POWER
% -------------------------------------------------------------------------
%    To use attenuation in the exponential field equations, we must convert 
%    dB/km to the linear attenuation coefficient alpha (m^-1).
%    
%    alpha [Np/m] = attenuation[dB/km] / (10 * log10(e) * 1000) [Ref 3, 4]
alpha = attenuation/(10 * log10(exp(1)) * 1000); 


%    To convert mW to dBm: 
%    P [dBm] = 10 * log10 (P [mW] / 1 [mW])
P_input_mW = 10; % Optical launch power (mW)
P_input_dBm = 10 * log10(P_input_mW / 1); % Same in dBm (10 dBm)


% -------------------------------------------------------------------------
% 1.4  SPATIAL DISCRETISATION
% -------------------------------------------------------------------------
dz = 0.05;      % Spatial sampling step (m): must satisfy dz << d
z = 0:dz:L-dz;  % Discrete position vector along the fiber (m)
Nz = length(z); % Total number of spatial samples


% -------------------------------------------------------------------------
% 1.5  PULSE-LIMITED SPATIAL RESOLUTION
% -------------------------------------------------------------------------
pulse_width = 10e-9; % Probe pulse duration (s) (10 ns)
d = c * pulse_width / (2 * n_ave);  % Two-way spatial resolution (m) (~1 m)
M = round(d / dz); % Number of scattering segments (inhomogeneities) within 
                   % one pulse


% -------------------------------------------------------------------------
% 1.6  LASER & NOISE PARAMETERS
% -------------------------------------------------------------------------
%    Laser linewidth: a narrower linewidth (e.g., 1 kHz) would increase 
%    coherence; a wider one (> 1 MHz) would significantly increase phase 
%    noise and degrade the correlation peak. [Ref 5]
linewidth = 100e3; % Laser linewidth (100 kHz)

%    SNR = 20 dB is a typical lab operating point; 15 dB adds realistic 
%    noise while keeping the algorithm testable.
SNR_dB = 15;       % Receiver signal-to-noise ratio (dB)
sigma_n = 2e-6;    % Standard deviation of background index fluctuations


% -------------------------------------------------------------------------
% 1.7  FREQUENCY SWEEP PARAMETERS
% -------------------------------------------------------------------------
freq_range = 1000e6;    % Total optical frequency scan range (Hz) (1 GHz)
delta_f = 5e6;          % Frequency tuning step (Hz) (5 MHz)

% Absolute frequency vector centered on the carrier nu0
f = nu0 + (-freq_range/2 : delta_f : freq_range/2);  
Nf = length(f);         % Number of frequency steps in the sweep


% -------------------------------------------------------------------------
% 1.8  PERTURBATION GEOMETRY
% -------------------------------------------------------------------------
num_events  = 5;        % Number of distinct sensing events to simulate
pert_length = 1.5;      % Axial length of each perturbation zone (m)
spacing = 4;            % Gap between consecutive events (m)
sensing_zone = L - 30;  % Starting position of the sensing region (m)


% -------------------------------------------------------------------------
% 1.9  DETECTION / POST-PROCESSING THRESHOLDS
% -------------------------------------------------------------------------
std_mult = 3;   % Noise threshold multiplier (x std of calibration trace)
match_tolerance = 2.0; % Spatial radius used to match detected peaks to 
                       % theoretical event positions (m)


% -------------------------------------------------------------------------
% 1.10  PRINT CONFIGURATION SUMMARY
% -------------------------------------------------------------------------
fprintf(' ==== PHI-OTDR SIMULATION CONFIGURATION ==== \n');
fprintf([ ...
    'Fiber Length (L):                      %d m\n' ...
    'Spatial Resolution (d):                %.2f m\n' ...
    'Total Spatial Points (Nz)              %d\n' ...
    'Number of Scattering Segments (M):     %d\n' ...
    'Sampling Interval (dz):                %.3f m\n' ...                   
    'Attenuation:                           %.2f dB/km\n' ...
    'Attenuation Coefficient (α):           %.4e m^-1\n' ...
    'Average Refractive Index (n_ave):      %.3f\n' ...
    'Operating Wavelength (lambda_0):       %d nm\n'...
    'Central frequency (nu_0)               %.4e Hz\n'...
    'Pulse Width:                           %d ns\n' ...
    'Input Power:                           %d mW (%.1f dBm)\n' ...
    'Sweep Range:                           %d MHz\n' ...
    'Frequency step (Δf):                   %d MHz\n' ...
    'Number of frequencies (Nf):            %d\n' ...
    'Laser Linewidth (Δν):                  %.1f kHz\n' ...
    'System SNR:                            %d dB\n' ...
    'Index standard deviation (σ_n):        %.1e\n' ...
    'Number of Events:                      %d\n'...
    'Width of each event:                   %.2f m\n' ...
    'Separation between events:             %.2f m\n' ...
    'Sensing Zone Start:                    %d m\n' ...
    'Spatial Validation Radius              %.1f m\n'], ...
    L, d, Nz, M, dz, attenuation, alpha, n_ave, lambda0*1e9, nu0, ...
    pulse_width*1e9, P_input_mW, P_input_dBm, freq_range/1e6, delta_f/1e6, ...
    Nf, linewidth/1e3, SNR_dB, sigma_n, num_events, pert_length, spacing, ...
    sensing_zone, match_tolerance);

%% =======================================================================
% 2. STOCHASTIC RAYLEIGH SCATTERING PROFILE
% ========================================================================
%   IN : sigma_n, n_ave, Nz, dz
%   OUT: r(1×Nz) :: local Fresnel reflection amplitude at each cell boundary
%        n(1×Nz) :: baseline refractive index profile
%
%   Each dz-step is treated as an interface between media whose index
%   deviates from n_ave by a small zero-mean Gaussian random variable.
%   The Fresnel formula (normal incidence) gives the reflection amplitude
%   at each boundary:  r(i) = (n(i) − n(i+1)) / (n(i) + n(i+1))
%
%   Modeling the fiber as a 1D waveguide with random inhomogeneities. 
%   According to the paper, Rayleigh scattering is simulated by small 
%   fluctuations in the refractive index along the fiber core.


% ------------------------------------------------------------------------
% 2.1  RANDOM INDEX FLUCTUATIONS (GUASSIAN)
% ------------------------------------------------------------------------
delta_n = sigma_n * randn(1, Nz);   % Zero-mean index perturbations
n = n_ave + delta_n;                % Composite refractive-index profile n(z)


% ------------------------------------------------------------------------
% 2.1 FRESNEL REFLECTION COEFFICIENT
% ------------------------------------------------------------------------
%   Calculation of the local reflection coefficient (r) at each interface.
%   The model treats each 'dz' step as a discrete boundary between media 
%   with slightly different refractive indices.
r = zeros(1, Nz);
for i = 1:Nz-1
    % Fresnel formula for normal incidence between two adjacent cells
    % This represents the local backscattering amplitude at position z
    r(i) = (n(i) - n(i+1)) / (n(i) + n(i+1));
end

%% =======================================================================
% 3. ENVIRONMENTAL PERTURBATION MODEL (STRAIN / TEMPERATURE)
% ========================================================================
%   IN: sensing_zone, num_events, pert_length, spacing, magnitudes, gamma, eta
%   OUT: delta_n_pert(1×Nz) :: index perturbation profile
%        n_pert(1×Nz) :: total index profile (baseline + perturbation)
%        theor_* arrays :: ground-truth labels for the performance report
%
%   Each event is a rectangular step change in refractive index superimposed
%   on n(z). Temperature variation is recovered from delta_n via:
%       delta_T = delta_n / (gamma + eta * n_ave)
%
%   Localized environmental perturbations (such as physical strain or 
%   temperature changes) are modeled as localized modulations applied to the 
%   fiber's refractive index profile.
%   Temperature variations shift the phase response via the thermo-optic 
%   coefficient (dn/dT). (For Silica glass: dn/dT is approximately 
%   1.1e-5 K^-1)

% Unique index change magnitudes per event 
%   positive = heating / tension,
%   negative = cooling / compression
magnitudes = [0.38e-7; 6.84e-7; 2.28e-7; -9.81e-7; -1.14e-7;];

% Validate that all events fit within the fiber
first_event = sensing_zone; 
last_event  = sensing_zone + (num_events-1)*(pert_length + spacing) ...
    + pert_length;
if last_event >= L
    error(['Perturbation events exceed fiber length. Reduce num_events or ' ...
        'sensing_zone.']);
end


delta_n_pert = zeros(1, Nz); % Initialise perturbation overlay to zero

% Pre-allocate ground-truth arrays for later comparison with detections
theor_starts  = zeros(1, num_events);
theor_ends    = zeros(1, num_events);
theor_shifts  = zeros(1, num_events); % Expected frequency shift (Hz)
theor_delta_n = zeros(1, num_events);
theor_delta_T = zeros(1, num_events);

fprintf('\n%s\n', repmat('=', 1, 69));
fprintf('                     PERTURBATION EVENTS - PHI-OTDR\n');
fprintf('%s\n', repmat('=', 1, 69));
fprintf('%-8s | %-16s | %-10s | %-9s | %-10s\n', 'Event', 'Location (m)', ...
    'Shift (MHz)','Delta_n', 'Delta_T (K)');
fprintf('%s\n', repmat('-', 1, 69));

for i = 1:num_events
    % Spatial boundaries of event i
    start_pos = 2 + sensing_zone + (i - 1) * (pert_length + spacing);
    end_pos = start_pos + pert_length;
    mag = magnitudes(i);

    % Spectral shift expected from this index change:
    %   Delta_nu = nu0 * delta_n / n_ave
    shift_MHz = (nu0 * mag / n_ave) / 1e6;

    % Map continuous positions to discrete grid indices (clamped to valid 
    % range)
    start_idx = max(1, round(start_pos / dz));
    end_idx = min(Nz, round(end_pos / dz));

    % Inject rectangular perturbation into the index overlay
    delta_n_pert(start_idx:end_idx) = mag;

    % Convert delta_n to an equivalent temperature change
    delta_T = temp_variation(mag, gamma, eta, n_ave);

    % Store ground-truth values
    theor_starts(i)  = start_pos;
    theor_ends(i)    = end_pos;
    theor_shifts(i)  = shift_MHz;
    theor_delta_n(i) = mag;
    theor_delta_T(i) = delta_T;

    loc_str = sprintf('[%.2f; %.2f]', start_pos, end_pos);
    fprintf('Event %-2d | %-16s | %+-11.2f | %+-9.2e | %+-10.4f\n', ...
        i, loc_str, shift_MHz, mag, delta_T);

end

% Superimpose perturbation on baseline index profile
n_pert = n + delta_n_pert;
fprintf('%s\n\n', repmat('-', 1, 69));

%% ========================================================================
% 4. PROBE SIGNAL & FREQUENCY SWEEP
% ========================================================================
%   IN : f, Nf, Nz, M, n_ave, dz, c
%   OUT: E_ref, E_sig       — noisy backscattered fields  (Nf × Nz-M+1)
%        E_ref_id, E_sig_id — noise-free reference fields (for benchmarking)
%        E_ref_raw_all      — pre-noise fields (for calibration noise floor)
%        z_valid            — trimmed spatial axis accounting for pulse width
%        t_laser            — round-trip time vector
%
%   In this stage, we simulate a frequency-swept probe signal to recover 
%   the Rayleigh Backscatter (RB) spectra, as detailed in the static 
%   measurement section of the paper.


% Trim the spatial axis: the sliding-window convolution loses M-1 samples
% at the far end, so valid positions run from z(1) to z(Nz-M+1).
z_valid = z(1 : (Nz - M + 1)); 

% Round-trip time of flight to each valid position
t_laser = (0:Nz-M) * (2 * n_ave * dz / c); 

% Pre-allocate backscattered field matrices
% - Rows represent frequency components; 
% - Columns represent spatial positions (traces)
% The length is Nz-M+1 because the pulse integration window 'M' reduces the 
% valid range.
E_ref = zeros(Nf, Nz - M + 1); % Noisy reference state
E_sig  = zeros(Nf, Nz - M + 1); % Noisy perturbed state
E_ref_id = zeros(Nf, Nz - M + 1); % Ideal (noise-free) reference
E_sig_id = zeros(Nf, Nz - M + 1); % Ideal (noise-free) perturbed
E_ref_raw_all = zeros(Nf, Nz - M + 1); % Raw reference before receiver noise
E_sig_raw_all = zeros(Nf, Nz - M + 1); % Raw perturbed before receiver noise


%% ------------------------------------------------------------------------
% 4.1  PULSE SHAPE SELECTION
% ------------------------------------------------------------------------
%   The pulse shape defines the convolution window 'window' (length M).
%   Two physically motivated models are offered:
%     1) RC-filter model :: first-order EOM bandwidth limit (fc = 100 MHz)
%     2) Super-Gaussian :: flat-top approximation of order N

fprintf('\nPulse shaping selection:     1) RC Filter Model      2) Super-Gaussian Model\n');
shape_choice = input('  > Option ');

if shape_choice == 1
    % ---------------------------------------------------------------------
    % 4.1.a) RC FILTER MODEL
    % ---------------------------------------------------------------------
    %   This section models the Electro-Optic Modulator (EOM) response as a 
    %   first-order RC low-pass filter. This simulates the hardware 
    %   constraints where the electronic driver cannot switch 
    %   instantaneously, resulting in finite rise and fall times
    % 
    %   Rise times follow V(t) = A(1−e^(−t/RC))
    %   Fall times follow V(t) = V(tau)·e^(−(t−tau)/RC)
        
    fc = 100e6;               % Filter cut-off frequency (Hz)
    RC = 1 / (2 * pi * fc);   % Time constant (s)
    tau_pulse = pulse_width;  % Ideal rectangular pulse duration (s)
    A = 1;                    % Normalized peak amplitude
    
    %   Time vector: Covers both the charging (rise) and discharging (fall) 
    %   phases. Although the physical pulse duration is tau_pulse, the RC
    %   tail extends beyond it. We generate 2*M points to visualize the 
    %   full waveform. 
    t_rc = linspace(0, 2 * pulse_width, 2 * M);
    
    % Ideal rectangular input for comparison
    rect_pulse = ones(1, M);
    
    % Initialize the full waveform vector
    window_full = zeros(1, 2 * M);
    
    % Apply the step-response equations for an RC circuit
    for k = 1:2*M
        tk = t_rc(k);
        if tk <= 0
            window_full(k) = 0;
        elseif tk < tau_pulse
            % Rising Edge (Charging Phase): 0 < t < tau
            % Equation: V(t) = A * (1 - exp(-t/RC))
            window_full(k) = A * (1 - exp(-tk / RC));
        else
            % Falling Edge (Discharging Phase): t > τ
            % The decay starts from the amplitude reached at the end of the pulse.
            % Equation: V(t) = V(tau) * exp(-(t - tau)/RC)
            window_full(k) = A * (1 - exp(-tau_pulse / RC)) * exp(-(tk - tau_pulse) / RC);
        end
    end
    
    % The 'window' used for the backscatter convolution represents the pulse 
    % within the sampling window M. 
    window = window_full(1:M);
    window = window / max(window); % Normalise peak to 1

    % ---------------------------------------------------------------------
    %                               VISUAL CHECK
    % ---------------------------------------------------------------------
    figure(1); set(gcf, 'Name', 'Pulse Shape');
    
    % Time vector for display purposes only (covering twice the pulse width)
    t_rc_ns = linspace(0, 2 * pulse_width * 1e9, 2 * M); % Time in ns
    t_window_ns = linspace(0, pulse_width * 1e9, M); % Time in ns for the 
                                                     % active window
    
    % Shows the full charging/discharging cycle of the RC filter
    plot(t_rc_ns, window_full, 'Color', [0.1 0.2 0.8], 'LineWidth', 2, ...
     'DisplayName', 'RC Pulse Shapping');
    hold on;
    plot(t_window_ns, ones(1,M), 'r--', 'LineWidth', 1.5, 'DisplayName', ...
        'Ideal Rectangular Pulse');
    hold off;

    xlabel('Time (ns)'); ylabel('Amplitude');
    title('RC Pulse Shapping (fc = 100 MHz)');
    legend('Location', 'northeast'); ylim([0 1.2]); grid on;
    
elseif shape_choice == 2
    % ---------------------------------------------------------------------
    % 4.1.b) SUPER-GAUSSIAN MODEL
    % ---------------------------------------------------------------------
    %   The Super-Gaussian function is used to model pulses that have a 
    %   smoother transition than a hard rectangular pulse but can maintain 
    %   a "flat-top" characteristic depending on the order (N).
    % 
    %   N = 1 :: standard Gaussian
    %   N >= 3 :: near-rectangular with realistic roll-off
    
    order_N = 3; % Super-Gaussian order
    
    % Symmetric time vector centred at zero
    t_sg = linspace(-pulse_width, pulse_width, M);
   
    % Generate the Super-Gaussian window using the sgauss function.
    % Parameters: (time, FWHM width, Energy, Chirp, Order)
    window = sgauss(t_sg, pulse_width, 1, 0, order_N);
    window = window / max(window); % Normalise peak to 1
    
    % Reference: Create a standard rectangular pulse for visual comparison.
    rect_pulse = ones(1, M); 

    % ---- VISUAL VERIFICATION (SUPER-GAUSSIAN SHAPE) ----
    figure(1); set(gcf, 'Name', 'Pulse Shape');
   
    plot(t_sg * 1e9, window, ...
        'Color', [0.2 0.2 0.8], 'LineWidth', 2, ...
        'DisplayName', ['Super-Gaussian (N = ' num2str(order_N) ')']);
    hold on;
    plot(t_sg * 1e9, rect_pulse, 'r--', 'LineWidth', 1.5, ...
         'DisplayName', 'Ideal Rectangular Pulse');
    hold off;
    title('Super Gaussian Model Pulse Shapping');
    xlabel('Relative Time (ns)'); ylabel('Normalized Amplitude');
    legend('Location', 'northeast'); ylim([0 1.2]); grid on;
end

%% ========================================================================
% 5. COHERENT BACKSCATTER INTEGRATION & NOISE MODELLING
% ========================================================================
%   IN : r, n, n_pert, f, Nf, c, dz, window, t_laser, alpha, z_valid,
%        P_input_dBm, linewidth, SNR_dB
%   OUT: E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_all, E_sig_raw_all
%
%   For each optical frequency, the backscattered electric field is built
%   by convolving the complex scattering profile 
%                             r(z) * exp(j * 2 * phi(z)) 
%   with the pulse window. 
%   'conv(..., ''valid'')' implements the sliding coherent summation over 
%   M scatterers — physically equivalent to the detector integrating all 
%   back-reflections that arrive within one pulse width.
%
%   Noise chain:
%     1. Laser phase noise :: Wiener process (random walk), scales with tau
%     2. Fiber loss :: Beer-Lambert field decay: exp(−alpha * z)
%     3. Receiver noise :: AWGN (thermal + shot) added via awgn()


for f_idx = 1:Nf
    current_nu = f(f_idx); % Current sweep frequency (Hz)
    shift = current_nu - nu0; % Offset from the carrier (Hz)
    
    % -------------------------------------------------------------------
    % 5.1 PROPAGATION PHASE (CUMULATIVE, ONE-WAY)
    % -------------------------------------------------------------------
    %   phi(z) = cumsum(beta(z))*dz  where beta(z) = 2*pi*n(z)*nu/c
    %   The term 'r .* exp(1j * 2 * phi)' represents the local 
    %   backscattered light from each segment, where '2*phi' accounts for 
    %   the round-trip path.
    
    % Reference state: Cumulative phase with baseline refractive index n
    beta_ref = 2 * pi * n * current_nu /c;
    phi_ref = cumsum(beta_ref) * dz; % Round-trip baseline phase vector

    % Perturbed state: Cumulative phase with modified index n_pert
    % This captures the phase shift induced by the sensing event (strain/temp).
    beta_sig = 2 * pi * n_pert * current_nu / c;
    phi_sig = cumsum(beta_sig) * dz; % Round-trip perturbed phase vector

    % -------------------------------------------------------------------
    % 5.2 COHERENT BACKSCATTER VIA SLIDING-WINDOW CONVOLUTION
    % -------------------------------------------------------------------
    %                           WHY CONVOLUTION ('conv')? 
    %   Physically, the detector captures the coherent phasor sum of all 
    %   M reflectors within the pulse width at any given time delay. 
    %   Mathematically, this is equivalent to a sliding window integration.
    %   Using 'conv' with a rectangular 'window' replaces a nested spatial 
    %   loop, significantly optimizing the simulation while maintaining 
    %   exact physical consistency with the 1D waveguide model.
    %
    %   Each sample of 'E_ref_conv' is the phasor sum of M reflectors
    %   illuminated by the pulse at that round-trip delay.
    E_ref_conv = conv(r .* exp(1j * 2 * phi_ref), window, 'valid');
    E_sig_conv = conv(r .* exp(1j * 2 * phi_sig), window, 'valid');

    % -------------------------------------------------------------------
    % 5.3 LASER PHASE NOISE (TRANSMITER IMPAIRMENT)
    % -------------------------------------------------------------------
    %   lasercw() returns a complex envelope whose phase executes a Wiener
    %   process with diffusion rate proportional to the laser linewidth.
    E_laser = lasercw(t_laser, P_input_dBm, 0, linewidth, shift);

    % -------------------------------------------------------------------
    % 5.4 ROUND-TRIP FIBER ATTENUATION (BEER_LAMBERT)
    % -------------------------------------------------------------------
    %   As the light travels to distance 'z' and back to the detector 
    %   (round-trip), the electric field amplitude decays exponentially.
    %   We use 'exp(-alpha * z)' because the signal accumulates loss over 
    %   the total path (2*z). Since alpha is defined for power, the field 
    %   decay over distance '2z' is exp(-(alpha/2) * 2z) = exp(-alpha * z).
    loss_factor = exp(-alpha * z_valid);

    % Combine scattering profile, phase noise, and loss
    E_ref_raw = E_ref_conv .* E_laser .* loss_factor;
    E_sig_raw = E_sig_conv .* E_laser .* loss_factor;

    % -------------------------------------------------------------------
    % 5.5 ADDITIVE RECEIVER NOISE (THERMAL + SHOT)
    % -------------------------------------------------------------------
    %   Simulates the electronic noise floor (Thermal) and photon counting 
    %   noise (Shot). The 'measured' flag ensures the noise power is scaled 
    %   relative to the signal power.
    E_ref(f_idx, :) = awgn(E_ref_raw, SNR_dB, 'measured');
    E_sig(f_idx, :) = awgn(E_sig_raw, SNR_dB, 'measured');

    % Store noise-free and pre-noise versions for post-processing analysis
    E_ref_id(f_idx, :) = E_ref_conv;
    E_sig_id(f_idx, :) = E_sig_conv;
    E_ref_raw_all(f_idx, :) = E_ref_raw;
    E_sig_raw_all(f_idx, :) = E_sig_raw;

end

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


lags_freq = (-(Nf-1):(Nf-1)) * delta_f; % Lag axis in Hz
% Pre-allocate matrix for 3D correlation visualization
corr_map = zeros(length(lags_freq), Nz-M+1);
% Vector to store the estimated frequency shift per position
freq_shift = zeros(1, Nz-M+1);               
freq_shift_calib = zeros(1, Nz - M +1);

for k = 1:Nz-M+1
    % Intensity spectra (power spectral density approximation)
    I_ref = abs(E_ref(:,k)).^2;
    I_sig = abs(E_sig(:,k)).^2;
    
    % Normalised cross-correlation of zero-mean spectra
    [cv, lags] = xcorr(I_sig - mean(I_sig), I_ref - mean(I_ref), 'coeff');
    corr_map(:, k) = cv;
    
    % Peak lag -> spectral shift at position k
    [~, max_idx] = max(cv);
    freq_shift(k) = lags(max_idx) * delta_f; 

    % Calibration path: correlate two independent realisations of the
    % reference state to estimate the noise-induced shift variance

    I_calib = abs(awgn(E_ref_raw_all(:,k), SNR_dB, 'measured')).^2;
    [cv_c, lags_c] = xcorr(I_calib - mean(I_calib), I_ref - mean(I_ref), 'coeff');
    [~, max_idx_c] = max(cv_c);
    freq_shift_calib(k) = lags_c(max_idx_c) * delta_f;
end
%% ========================================================================
% 7. PEAK DETECTION WITH NOISE SUPPRESSION
% ========================================================================
%   IN : freq_shift, freq_shift_calib, z_valid, M, dz, pert_length, std_mult
%   OUT: locs_pos, pks_pos — positions and amplitudes of positive peaks
%        locs_neg, pks_neg — positions and amplitudes of negative peaks
%        all_locs, all_pks — merged, sorted event list
%        smooth_freq_shift — low-pass smoothed shift profile
%        threshold         — adaptive noise floor

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


% -------------------------------------------------------------------
% 7.2 ADAPTIVE NOISE THRESHOLD
% -------------------------------------------------------------------
%   The standard deviation of the calibration trace reflects the noise
%   floor; events must exceed std_mult*sigma(freq_shift_calib) to be 
%   declared.
threshold = std_mult * std(freq_shift_calib); 


% -------------------------------------------------------------------
% 7.3 SEPARATED POSITVE/NEGATIVE PEAK SEARCH
% -------------------------------------------------------------------
%   Events with opposite-sign delta_n produce opposite-sign shifts, so
%   positive and negative peaks are detected independently.
min_peak_dist = round(pert_length/dz); % Minimum inter-peak grid spacing

% Common findpeaks settings:
%   MinPeakHeight :: rejects all variations below the threshold floor
%   MinPeakDistance :: ensures detected peaks are separated by at least one 
%                      event width
%   MinPeakProminence :: measures peak height relative to the surrounding 
%                        baseline, filtering out ripple artifacts located 
%                        on long signal slopes
peak_opts = { ...
    'MinPeakHeight',      threshold, ...
    'MinPeakDistance',    min_peak_dist * dz, ...
    'MinPeakProminence',  threshold * 0.3 };

[pks_pos, locs_pos] = findpeaks( smooth_freq_shift, z_valid, peak_opts{:});
[pks_neg, locs_neg] = findpeaks(-smooth_freq_shift, z_valid, peak_opts{:});
pks_neg = -pks_neg; % Restore true (negative) amplitudes

% -------------------------------------------------------------------
% 7.4 MERGE AND STOR ALL DETECTED EVENTS BY POSITION
% -------------------------------------------------------------------
all_locs = [locs_pos, locs_neg];
all_pks  = [pks_pos,  pks_neg];
[all_locs, sort_idx] = sort(all_locs);
all_pks = all_pks(sort_idx);

%% ========================================================================
% 8. PERFORMANCE REPORT
% ========================================================================
%   Calls the external 'report' utility which prints a confusion matrix,
%   per-event detection status, and figures-of-merit (precision, recall, F1).
report(theor_starts, theor_ends, theor_shifts, theor_delta_n, ...
    theor_delta_T, all_locs, all_pks, z_valid, smooth_freq_shift, ...
    threshold, gamma, eta, n_ave, nu0, pert_length, match_tolerance);

%% ========================================================================
% 9. DATA VISUALISATION
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
plot(z_valid, 10*log10(P_ref_ideal + eps), 'b', 'DisplayName', ...
    'Ideal (No Loss)');
hold on;
plot(z_valid, 10*log10(P_ref_loss + eps), 'r', 'DisplayName', ...
    ['Fiber Loss (', num2str(attenuation), ' dB/km) + Noise']);
hold off;
title('Backscattered Intensity (Logarithmic Scale)'); 
ylabel('Power (dBm)'); xlabel('Distance (m)'); grid on;  legend('Location', 'southeast');

%%
% --- FIGURE 7: Peak Detection Diagnostics ---
% Visualizes raw versus smoothed frequency shifts alongside thresholds and 
% positive/negative peak markers.
figure(7);
set(gcf, 'Name', 'Event Detection');
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

%%

fprintf('\n--- Simulation successfully completed! ---\n');

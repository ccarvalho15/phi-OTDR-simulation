function conf = get_configuration()
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
    conf.c = 3e8;                    % Speed of light in vacuum (m/s)
    conf.lambda0 = 1550e-9;          % Operating wavelength (m)
    conf.nu0 = conf.c/conf.lambda0;            % Central optical frequency (Hz)
    
    
    % -------------------------------------------------------------------------
    % 1.2  FIBER MEDIUM PROPERTIES
    % -------------------------------------------------------------------------
    conf.L = 240;           % Total fiber length used in the lab (m)
    conf.n_ave = 1.456;     % Average refractive index of the silica core
    conf.gamma = 9.1e-6;    % Thermo-optic coefficient  dn/dT (K^-1)
    conf.eta = 0.5e-6;      % Thermal expansion coefficient (K^-1)
    
    %    Attenuation is typically provided in dB/km (logarithmic scale). Standard 
    %    single-mode fiber (SMF-28) at 1550nm has approx. 0.19–0.2 dB/km range. 
    %    [Refs 3,4]
    conf.attenuation = 0.2; % Fiber attenuation (dB/km)
    
    
    % -------------------------------------------------------------------------
    % 1.3  LOSS CONVERSION & INPUT POWER
    % -------------------------------------------------------------------------
    %    To use attenuation in the exponential field equations, we must convert 
    %    dB/km to the linear attenuation coefficient alpha (m^-1).
    %    
    %    alpha [Np/m] = attenuation[dB/km] / (10 * log10(e) * 1000) [Ref 3, 4]
    conf.alpha = conf.attenuation/(10 * log10(exp(1)) * 1000); 
    
    
    %    To convert mW to dBm: 
    %    P [dBm] = 10 * log10 (P [mW] / 1 [mW])
    conf.P_input_mW = 10; % Optical launch power (mW)
    conf.P_input_dBm = 10 * log10(conf.P_input_mW / 1); % Same in dBm (10 dBm)
    
    
    % -------------------------------------------------------------------------
    % 1.4  SPATIAL DISCRETISATION
    % -------------------------------------------------------------------------
    conf.dz = 0.05;      % Spatial sampling step (m): must satisfy dz << d
    conf.z = 0:conf.dz:conf.L-conf.dz;  % Discrete position vector along the fiber (m)
    conf.Nz = length(conf.z); % Total number of spatial samples
    
    
    % -------------------------------------------------------------------------
    % 1.5  PULSE-LIMITED SPATIAL RESOLUTION
    % -------------------------------------------------------------------------
    conf.pulse_width = 10e-9; % Probe pulse duration (s) (10 ns)
    conf.d = conf.c * conf.pulse_width / (2 * conf.n_ave);  % Two-way spatial resolution (m) (~1 m)
    conf.M = round(conf.d / conf.dz); % Number of scattering segments (inhomogeneities) within 
                       % one pulse
    
    
    % -------------------------------------------------------------------------
    % 1.6  LASER & NOISE PARAMETERS
    % -------------------------------------------------------------------------
    %    Laser linewidth: a narrower linewidth (e.g., 1 kHz) would increase 
    %    coherence; a wider one (> 1 MHz) would significantly increase phase 
    %    noise and degrade the correlation peak. [Ref 5]
    conf.linewidth = 100e3; % Laser linewidth (100 kHz)
    
    %    SNR = 20 dB is a typical lab operating point; 15 dB adds realistic 
    %    noise while keeping the algorithm testable.
    conf.SNR_dB = 15;       % Receiver signal-to-noise ratio (dB)
    conf.sigma_n = 2e-6;    % Standard deviation of background index fluctuations
    
    
    % -------------------------------------------------------------------------
    % 1.7  FREQUENCY SWEEP PARAMETERS
    % -------------------------------------------------------------------------
    conf.freq_range = 1000e6;    % Total optical frequency scan range (Hz) (1 GHz)
    conf.delta_f = 5e6;          % Frequency tuning step (Hz) (5 MHz)
    
    % Absolute frequency vector centered on the carrier nu0
    conf.f = conf.nu0 + (-conf.freq_range/2 : conf.delta_f : conf.freq_range/2);  
    conf.Nf = length(conf.f);         % Number of frequency steps in the sweep
    
    
    % -------------------------------------------------------------------------
    % 1.8  PERTURBATION GEOMETRY
    % -------------------------------------------------------------------------
    conf.num_events  = 5;        % Number of distinct sensing events to simulate
    conf.pert_length = 1.5;      % Axial length of each perturbation zone (m)
    conf.spacing = 4;            % Gap between consecutive events (m)
    conf.sensing_zone = conf.L - 30;  % Starting position of the sensing region (m)
    
    
    % -------------------------------------------------------------------------
    % 1.9  DETECTION / POST-PROCESSING THRESHOLDS
    % -------------------------------------------------------------------------
    conf.std_mult = 3;   % Noise threshold multiplier (x std of calibration trace)
    conf.match_tolerance = 2.0; % Spatial radius used to match detected peaks to 
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
        conf.L, conf.d, conf.Nz, conf.M, conf.dz, conf.attenuation, ...
        conf.alpha, conf.n_ave, conf.lambda0*1e9, conf.nu0, ...
        conf.pulse_width*1e9, conf.P_input_mW, conf.P_input_dBm, ...
        conf.freq_range/1e6, conf.delta_f/1e6, conf.Nf, conf.linewidth/1e3, ...
        conf.SNR_dB, conf.sigma_n, conf.num_events, conf.pert_length, ...
        conf.spacing, conf.sensing_zone, conf.match_tolerance);
end


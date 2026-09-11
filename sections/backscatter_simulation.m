function [z_valid, t_laser, E_ref, E_sig, E_ref_id, E_sig_id, E_ref_raw_all, E_ref_raw_calib_all] = ...
    backscatter_simulation(conf, shape_choice, r, n, n_pert)

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

z           = conf.z;
Nz          = conf.Nz;
M           = conf.M;
n_ave       = conf.n_ave;
dz          = conf.dz;
c           = conf.c;
Nf          = conf.Nf;
f           = conf.f;
nu0         = conf.nu0;
alpha       = conf.alpha;
pulse_width = conf.pulse_width;
P_input_dBm = conf.P_input_dBm;
linewidth   = conf.linewidth;
SNR_dB      = conf.SNR_dB;

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
E_ref_raw_calib_all = zeros(Nf, Nz - M + 1);


%% ------------------------------------------------------------------------
% 4.1  PULSE SHAPE SELECTION
% ------------------------------------------------------------------------
%   The pulse shape defines the convolution window 'window' (length M).
%   Two physically motivated models are offered:
%     1) RC-filter model :: first-order EOM bandwidth limit (fc = 100 MHz)
%     2) Super-Gaussian :: flat-top approximation of order N
%     3) Rectangular Pulse :: Ideal pulse with instantaneous transitions

if shape_choice == 1
    % ---------------------------------------------------------------------
    % 4.1.a) RC FILTER MODEL
    % ---------------------------------------------------------------------
    fc = 100e6;               % Filter cut-off frequency (Hz)
    RC = 1 / (2 * pi * fc);   % Time constant (s)
    tau_pulse = pulse_width;  % Ideal rectangular pulse duration (s)
    A = 1;                    % Normalized peak amplitude
    
    t_rc = linspace(0, 2 * pulse_width, 2 * M);
    window_full = zeros(1, 2 * M);
    
    for k = 1:2*M
        tk = t_rc(k);
        if tk <= 0
            window_full(k) = 0;
        elseif tk < tau_pulse
            window_full(k) = A * (1 - exp(-tk / RC));
        else
            window_full(k) = A * (1 - exp(-tau_pulse / RC)) * exp(-(tk - tau_pulse) / RC);
        end
    end
    
    window = window_full(1:M);
    window = window / max(window); % Normalise peak to 1

    % ---------------------------------------------------------------------
    %                               VISUAL CHECK
    % ---------------------------------------------------------------------
    figure(1); set(gcf, 'Name', 'Pulse Shape');
    
    % Centrar vetor de tempo RC para visualização no intervalo -20 a 20 ns
    t_rc_ns = linspace(0, 2 * pulse_width * 1e9, 2 * M); 
    
    plot(t_rc_ns, window_full, 'Color', [0.1 0.2 0.8], 'LineWidth', 3, ...
         'DisplayName', 'RC Pulse Shaping (f_c = 100 MHz)');
    hold on;
    
    % Ideal rectangular pulse
    t_rect_ns = [0, 0, pulse_width*1e9, pulse_width*1e9];
    v_rect    = [0, 1, 1, 0];
    plot(t_rect_ns, v_rect, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Ideal Rectangular Pulse');
    hold off;

    % Ajuste de fontes e eixos
    xlim([-5, 25]); ylim([0 1.2]); grid on;
    set(gca, 'FontSize', 20, 'LineWidth', 1.5); % Aumenta números dos eixos
    
    xlabel('Relative Time (ns)', 'FontSize', 20); 
    ylabel('Normalized Amplitude', 'FontSize', 20);
    title('RC Low-Pass Filtered Pulse Shaping', 'FontSize', 22, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 15);

elseif shape_choice == 2
    % ---------------------------------------------------------------------
    % 4.1.b) SUPER-GAUSSIAN MODEL
    % ---------------------------------------------------------------------
    order_N = 3; % Super-Gaussian order

    t_sg = linspace(0, pulse_width, M);
    t_center = pulse_width/2;

    window = sgauss(t_sg - t_center, pulse_width, 1, 0, order_N);
    window = window / max(window); % Normalise peak to 1

    % ---------------------------------------------------------------------
    %                               VISUAL CHECK
    % ---------------------------------------------------------------------
    figure(1); set(gcf, 'Name', 'Pulse Shape');

    t_plot = linspace(-pulse_width, 2*pulse_width, 2*M);
    window_plot = sgauss(t_plot - t_center, pulse_width, 1, 0, order_N);
    window_plot = window_plot / max(window_plot);

    plot(t_plot * 1e9, window_plot, 'Color', [0.2 0.2 0.8], 'LineWidth', 3, ...
         'DisplayName', ['Super-Gaussian (p = ' num2str(order_N) ')']);
    hold on;
    
    % Ideal rectangular pulse
    t_rect_ns = [0, 0, pulse_width*1e9, pulse_width*1e9];
    v_rect    = [0, 1, 1, 0];
    plot(t_rect_ns, v_rect, 'r--', 'LineWidth', 2.5, 'DisplayName', 'Ideal Rectangular Pulse');
    hold off;

    % Ajuste de fontes e eixos
    xlim([-15, 25]); ylim([0 1.2]); grid on;
    set(gca, 'FontSize', 20, 'LineWidth', 1.5); % Aumenta números dos eixos
    
    xlabel('Relative Time (ns)', 'FontSize', 20); 
    ylabel('Normalized Amplitude', 'FontSize', 20);
    title('Super Gaussian Pulse Shapping', 'FontSize', 22, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 15);

elseif shape_choice == 3
    % -------------------------------------------------------------------------
    % 4.1.c) IDEAL RECTANGULAR PULSE
    % -------------------------------------------------------------------------
    window = ones(1, M);
    window = window / max(window); % Normalise peak to 1

    % ---------------------------------------------------------------------
    %                               VISUAL CHECK
    % ---------------------------------------------------------------------
    figure(1); set(gcf, 'Name', 'Pulse Shape');
    
    t_rect_ns = [0, 0, pulse_width*1e9, pulse_width*1e9];
    v_rect    = [0, 1, 1, 0];
    
    plot(t_rect_ns, v_rect, 'r-', 'LineWidth', 3, 'DisplayName', 'Ideal Rectangular Pulse');

    % Ajuste de fontes e eixos
    xlim([-11, 21]); ylim([0 1.2]); grid on;
    set(gca, 'FontSize', 18, 'LineWidth', 1.5); % Aumenta números dos eixos
    
    xlabel('Relative Time (ns)', 'FontSize', 18); 
    ylabel('Normalized Amplitude', 'FontSize', 18);
    title('Ideal Rectangular Pulse', 'FontSize', 20, 'FontWeight', 'bold');
    % legend('Location', 'northeast', 'FontSize', 16);
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

    % -------------------------------------------------------------------
    % 5.6 CALIBRATION TRACE — INDEPENDENT PHASE-NOISE REALIZATION
    % -------------------------------------------------------------------
    %   Generated AFTER the primary signal (E_ref, E_sig) is finalized,
    %   so this extra random draw does not shift the noise stream used
    %   for the primary reference/perturbed fields. Only the calibration
    %   trace — and therefore the adaptive threshold — is affected.
    E_laser_calib = lasercw(t_laser, P_input_dBm, 0, linewidth, shift);
    E_ref_raw_calib = E_ref_conv .* E_laser_calib .* loss_factor;
    E_ref_raw_calib_all(f_idx, :) = E_ref_raw_calib;
end

end
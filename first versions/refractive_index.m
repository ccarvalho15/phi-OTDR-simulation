%% REFRACTIVE INDEX CHANGE SIMULATION
close all; clear; clc;

%% ---------------- PARAMETERS ----------------
c = 3e8;                    % Speed of light (m/s)
lambda0 = 1550e-9;          % Central wavelength (m)
nu0 = c/lambda0;            % Central optical frequency (Hz)
n_ave = 1.456;              % Average refractive index of the silica fiber

L = input('> Total fiber length (m), L: ');                   
dz = input('> Spatial sampling interval (m), dz: ');
z = 0:dz:L-dz;              % Distance vector
Nz = length(z);

pulse_width = 10e-9;        % Pulse width (10 ns)
d = c*pulse_width/(2*n_ave);% Spatial resolution (~1 m)
M = round(d/dz);            % Number of points contained within one pulse width

if M < 10
    warning('At least 10 samples per pulse are required (M = %d).', M);
end

alpha = 0;                  % Fiber attenuation
gamma = 9.1e-6;             % Thermo-optic coefficient
eta = 0.5e-6;               % Thermal expansion coefficient

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

%% ------------ RANDOM TEMPERATURE PERTURBATIONS ------------
n_pert = n;
num_events = input('> Enter number of events to simulate: ');
width_idx = round(1 / dz); % FWHM of 1 meter converted to indices

if num_events > 1
    fixed_dist = input('> Enter the fixed distance between events (m): ');
    step_idx = round(fixed_dist / dz);
    % Ensure the last event does not exceed Nz-M
    max_start_idx = (Nz - M) - (num_events - 1) * step_idx;
    if max_start_idx <= M
        error('Fiber is too short for the requested number of events!');
    end
else
    fixed_dist = 0; step_idx = 0; max_start_idx = Nz-M;
end

first_pos_idx = randi([M, max_start_idx]);
event_indices = first_pos_idx + (0:num_events-1) * step_idx;

delta_T = zeros(1, Nz);
fprintf('\n--- GROUND TRUTH: Perturbation Locations ---\n');

for i = 1:num_events
    rand_dt = (rand - 0.5) * 2.0; % Random Delta T (-1 K to 1 K)

    pos_idx = event_indices(i); 

    % Define the influence range of the event (e.g., 1 meter)
    range_idx = pos_idx : min(pos_idx + width_idx, Nz);
    delta_T(range_idx) = rand_dt;

    fprintf('Event %d: Position %.2f m | dT: %.3f K\n', ...
                i, z(pos_idx), rand_dt);
end
fprintf('-------------------------------------------\n');

% Calculation of the new refractive index under perturbation
n_pert = n + (gamma + n_ave * eta) * delta_T;

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
for f_idx = 1:Nf
    current_nu = f(f_idx);

    % Reference propagation phase (no change)
    beta_ref = 2 * pi * n * current_nu /c;
    phi_ref = cumsum(beta_ref)*dz;      % Integrated phase

    % Signal propagation phase (with temperature change)
    beta_sig = 2 * pi * n_pert * current_nu / c;
    phi_sig = cumsum(beta_sig)*dz;

    for k = 1:Nz-M
        idx = k:(k+M-1);
        % Complex field summation for both reference and measurement traces
        E_ref(f_idx, k) = sum(r(idx) .* exp(1j * 2 * phi_ref(idx)));
        E_sig(f_idx, k) = sum(r(idx) .* exp(1j * 2 * phi_sig(idx)));
    end
end

%% ------------ CROSS-CORRELATION ------------
freq_shift = zeros(1, Nz-M);
lags_freq = (-(Nf-1):(Nf-1)) * delta_f;
corr_map = zeros(length(lags_freq), Nz-M);

for k = 1:Nz-M
    ref = abs(E_ref(:,k)).^2;
    sig = abs(E_sig(:,k)).^2;

    % Normalization for zero-mean cross-correlation
    ref_norm = ref - mean(ref); 
    sig_norm = sig - mean(sig);

    % Cross-correlation of zero-mean intensity spectra
    [cv, lags] = xcorr(sig_norm, ref_norm, 'coeff');
    corr_map(:, k) = cv;

    % Finding the peak index to determine the frequency shift
    [~, max_idx] = max(cv);
    freq_shift(k) = lags(max_idx) * delta_f; 
    
end

%% ------------ DETECTION RESULTS IN CONSOLE ---------------
fprintf('\n--- DETECTION RESULTS ---');
threshold = 5e6; % 5 MHz threshold to ignore noise
[pks, locs] = findpeaks(abs(freq_shift), 'MinPeakHeight', threshold, 'MinPeakDistance', M);

if isempty(locs)
    fprintf('\nNo significant events detected above %.1f MHz threshold.', threshold/1e6);
else
    for n = 1:length(locs)
        fprintf('\nEvent %d: Position: %.2f m | Shift: %.2f MHz | Corr: %.4f', n, ...
            z(locs(n)), freq_shift(locs(n))/1e6, max(corr_map(:, locs(n))) );
    end
    [max_val, max_idx] = max(abs(freq_shift));
end
fprintf('\n---------------------------------\n');


%% ------------ 3D VISUALIZATION ------------
fprintf('\nGenerating multi-view 3D visualization...');

% Automatic limits to focus on the perturbation area
if ~isempty(locs)
    % Apply a -10m and +10m buffer around detected events
    dist_min = min(z(locs)) - 10;
    dist_max = max(z(locs)) + 10;

    % Security check: ensure limits remain within the fiber length [0, L]
    dist_min = max(0, dist_min);
    dist_max = min(L, dist_max);
else
    % Fallback if no events are detected
    dist_min = 0;
    dist_max = L;
end

f_min = -freq_range/2e6; f_max = freq_range/2e6;

% Event description for the title
if num_events == 1
    event_type_str = 'Single Event';
else
    event_type_str = sprintf('%d Events (Spacing: %.1f m)', num_events, fixed_dist);
end

% Construct a clean figure name for the window title bar
fig_name = sprintf('Phi-OTDR | %s | L=%dm, dz=%.2fm', event_type_str, L, dz);

% Create the figure
figure('Name', fig_name, 'NumberTitle', 'off');

events = length(locs);
colors = lines(events);

[Z_mesh, F_mesh] = meshgrid(z(1:Nz-M), lags_freq / 1e6);

% Define plot views and titles
views = [35, 75; 0, 90; 135, 75; 90, 0];
titles = {'Perspective View', 'Top View (XY)', 'Rear View', 'Side View (Frequency Shift)'};

h_plots = [];
legend_labels = cell(1, events);

for i = 1:4
    subplot(2, 2, i);
    
    % Surface Plot of Correlation Map
    surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none');
    colormap('jet'); colorbar; hold on;

    % Annotate detected events
    for n = 1:events
        dist_det = z(locs(n));
        shift_det = freq_shift(locs(n)) / 1e6;

        % Vertical indicator line
        plot3([dist_det dist_det], [shift_det shift_det], [-0.5 1], '--', ...
            'Color', colors(n,:), 'LineWidth', 1.5);

        % Peak marker
        p = plot3(dist_det, shift_det, 1, 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors(n,:), 'MarkerSize', 8);

        h_plots(n) = p;
        legend_labels{n} = sprintf('Event %d @ %.1f m, %.1f MHz', n, dist_det, shift_det);
    end

    view(views(i,1), views(i,2));
    xlabel('Distance (m)'); ylabel('Frequency Shift (MHz)'); zlabel('Cross-correlation (a.u.)')
    title(titles{i});

    % Axis Limits
    xlim([dist_min dist_max]); % Distance
    ylim([f_min f_max]); % Frequency shift
    zlim([-0.5 1.1]); % Correlation value
    grid on
end

% Create one global legend for the entire figure
if ~isempty(h_plots)
    llgd = legend(h_plots, legend_labels, 'Location', 'best');
end
rotate3d on;

fprintf(' Done!\n');
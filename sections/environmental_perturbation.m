function [gt, n_pert] = environmental_perturbation(conf, n, magnitudes)

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
    
    sensing_zone = conf.sensing_zone;
    num_events   = conf.num_events;
    pert_length  = conf.pert_length;
    spacing      = conf.spacing;
    Nz           = conf.Nz;
    dz           = conf.dz;
    nu0          = conf.nu0;
    n_ave        = conf.n_ave;
    gamma        = conf.gamma;
    eta          = conf.eta;
    L            = conf.L;
   
    % Validate that all events fit within the fiber
    last_event  = sensing_zone + (num_events-1)*(pert_length + spacing) ...
        + pert_length;
    if last_event >= L
        error(['Perturbation events exceed fiber length. Reduce num_events or ' ...
            'sensing_zone.']);
    end
    
    
    delta_n_pert = zeros(1, Nz); % Initialise perturbation overlay to zero
    
    % Pre-allocate ground-truth arrays for later comparison with detections
    gt.theor_starts  = zeros(1, num_events);
    gt.theor_ends    = zeros(1, num_events);
    gt.theor_shifts  = zeros(1, num_events); % Expected frequency shift (Hz)
    gt.theor_delta_n = zeros(1, num_events);
    gt.theor_delta_T = zeros(1, num_events);
    
    % fprintf('\n%s\n', repmat('=', 1, 69));
    % fprintf('                     PERTURBATION EVENTS - PHI-OTDR\n');
    % fprintf('%s\n', repmat('=', 1, 69));
    % fprintf('%-8s | %-16s | %-10s | %-9s | %-10s\n', 'Event', 'Location (m)', ...
    %     'Shift (MHz)','Delta_n', 'Delta_T (K)');
    % fprintf('%s\n', repmat('-', 1, 69));
    
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
        gt.theor_starts(i)  = start_pos;
        gt.theor_ends(i)    = end_pos;
        gt.theor_shifts(i)  = shift_MHz;
        gt.theor_delta_n(i) = mag;
        gt.theor_delta_T(i) = delta_T;
    
        loc_str = sprintf('[%.2f; %.2f]', start_pos, end_pos);
        % fprintf('Event %-2d | %-16s | %+-11.2f | %+-9.2e | %+-10.4f\n', ...
        %     i, loc_str, shift_MHz, mag, delta_T);
    
    end
    
    % Superimpose perturbation on baseline index profile
    n_pert = n + delta_n_pert;
    % fprintf('%s\n\n', repmat('-', 1, 69));
end
function [r, n] = rayleigh_profile(conf)

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

    sigma_n = conf.sigma_n;
    n_ave   = conf.n_ave;
    Nz      = conf.Nz;

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

end





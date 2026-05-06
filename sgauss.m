function y = sgauss(t,Tfwhm,E,C,m)
% Super-Gaussian pulse of complex amplitude centered at zero time
% USAGE:
% y = sgauss(t);
% y = sgauss(t,Tfwhm);
% y = sgauss(t,Tfwhm,E);
% y = sgauss(t,Tfwhm,E,C,m);
% y = sgauss(t,Tfwhm,E,C,m);
% INPUT
% t         time 
% Tfwhm     full-width at half maximum of the pulse power (default = 1)
% E         pulse energy (default = 1)
% C         chirp parameter (default = 0 for unchirped pulse)
% m         pulse order (sharpness) (default = 1 for Gaussian shape)
% OUTPUT
% y         vector of the same size as t, representing pulse amplitude A(t)
% PULSE DESCRIPTION
% A(t) = sqrt(Pin)*exp(-1/2*(1 + 1i*C)*(t/T0).^(2*m))
% Pin is the peak input power 
% Instantaneous power = |A(t)|^2
% Instantaneous phase = -(C/2)*(t/T0)^(2*m)
% NOTE
%   The time parameter T0 is related to Tfwhm by
%   T0 = Tfwhm/(2*(log(2)^(1/(2*m)))); 
%   The 10%-90% risetime Tr is given by
%   Tr = T0*(log(10)^(1/2m)- (log(10/9)^(1/2m))
%   The pulse peak power Pin = E/((T0/m)*gamma(1/(2*m))), where gamma is
%   the gamma function

% EXAMPLE
% dt = 0.01; % time step
% t = -2:dt:2;
% E = 1;
% Tfwhm = 1;
% C = 1;  
% m = 2; 
% y = sgauss(t,Tfwhm,E,C,m);
% phi = angle(y); % phase angle 
% subplot(3,1,1);
% plot(t,abs(y).^2); % plot instantaneous power
% title('2nd order Super Gaussian pulse: E = 1, Tfwhm = 1, C = 1')
% xlabel('Normalized time'); ylabel('Power');
% subplot(3,1,2)
% plot(t,phi) % plot instantaneous phase (radians)
% xlabel('Normalized time'); ylabel('phase (rad)');
% chirp = 1/(2*pi)*gradient(unwrap(phi),t); % plot chirp (Hz)
% subplot(3,1,3)
% plot(t,chirp) % plot chirp
% xlabel('Normalized time'); ylabel('chirp = d\phi/dt');
% PulseEnergy = trapz(abs(y).^2)*dt  % verify that pulse energy is equal to E

% By Prof. Michael Connelly
% Dept. Electronic and Computer Engineering
% University of Limerick, IRELAND
% michael.connelly@ul.ie
% Version 2.0 (April 2017)

if (nargin<5)
   m = 1;
end
if (nargin<4)
    C = 0; 
end
if (nargin<3)
  E = 1;
end
if (nargin<2)
  Tfwhm = 1;
end

T0 = Tfwhm/(2*(log(2)^(1/(2*m)))); 
P = E/((T0/m)*gamma(1/(2*m)));  % peak pulse power
% E is equal to integral of |A(t)|^2 over all time
y = sqrt(P)*exp(-((1 + 1j*C)/2)*(t/T0).^(2*m));

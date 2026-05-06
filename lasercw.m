function Eout=lasercw(t,powerdBm,initPhase,linewidth,varargin)
%% LASERCW - Continuous Wave Laser
%   The linewidth of the generated CW signal is modeled using a Gaussian
%   white noise source with variance of 2pi*linewidth corresponding to the
%   laser FWHM. The states of polarization is not considered in this model.
%
%
% Input arguments:
%              t - Time vector (starting at zero with constant step)[s]
%       powerdBm - Average output power [dBm]
%      initPhase - Initial phase [rad]
%      linewidth - FWHM Linewidth (Hz)
%     varargin:
% deltaFrequency - Frequency [Hz]
%  noiseResetGen - Flag to reset the random number generator{0,1}
%     noiseState - State of random number generator[]
%
% Output arguments:
%       Eout - Output field
%
% Call-examples:
%        Eout=lasercw(t,powerdbm,initPhase)
%        Eout=lasercw(t,powerdbm,initPhase,linewidth,resetNoiseGen,noiseState)

% Reference:
%
%[1]Eric Alpman, Florent Munier, Thomas Eriksson, Arne Sevensson, and Herbert
%  Airath. Estimation of phase noise for QPSK modulation over AWGN channels.
%  http://www.ep.liu.se/ecp/008/posters/004/ecp00804p.pdf.
%
%[2]Xiaopei Chen "Ultra-Narrow Laser Linewidth Measurement",PhD
%   Thesis,Virginia Polytechnic Institute and State University,july 2006.

%/////////////////////////////////////////////////////////////%
%(c)Carmo Medeiros                                        %
% Routine developed for : Radio-Over-Fiber Toolbox            %
% Last Revision: 10-10-2012                                  %
%/////////////////////////////////////////////////////////////%

global PARA

% Sample rate
Ts=t(2)-t(1);
fs=1/Ts;

if nargin <4
    error('The input arguments are not enough')
elseif nargin>7
    error('The input arguments are too many')
end
if nargin>6
    noiseState = varargin{3};
end
if nargin>5
    noiseResetGen = varargin{2};
end
if nargin>4
  deltaFrequency = varargin{1};
end

if nargin <=5
    noiseResetGen=0;
    noiseState=[];
end
if nargin <=4
    noiseResetGen=0;
    noiseState=[];
    deltaFrequency=0;
end



%% add linewidth noise
if linewidth>0
    
    if noiseResetGen==1
        rng('shuffle'); % Seed based on the current time. randn returns different values each time you do this.
    end
    
    % white noise variance = 2*pi*LW
    
    %% 1ºMethod
    
    noiseVariance = 2*pi*linewidth/(fs/2);
    
    if ~isempty(noiseState)
        rng(noiseState)
    end
    
    noise = randn(size(t))*sqrt(noiseVariance);
    
    phaseNoise=cumtrapz(noise);
    phaseNoise(end)=phaseNoise(1);
    
    %% 2º Method
    %     phaseNoise(1)=0;
    %     for i=2:length(t)
    %         phaseNoise(i)=phaseNoise(i-1)+randn(1,1)*sqrt(2*pi*linewidth*Ts);
    %     end
    %     phaseNoise(end)=0;
    %===
    
    %% debug
    %     figure
    %     plot(t/1e-9,phaseNoise/pi)
    %     xlabel('Time[ns]'),ylabel('Normalized phase by \pi')
    
    FM = exp(1i*phaseNoise);
    
elseif linewidth==0
    FM = 1;
end

%% create signal
powermW = 10^(powerdBm/10);

if deltaFrequency == 0
    Eout = sqrt(powermW*1e-3).*FM*exp(1i*initPhase);%[V]
else
    
    Eout = sqrt(powermW*1e-3).*FM.*exp(1i*2*pi*deltaFrequency*t).*exp(1i*initPhase);%[V]
end


end
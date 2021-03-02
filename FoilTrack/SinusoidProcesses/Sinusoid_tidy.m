%% Sinusoid_tidy.
% Makes phase of sinewave to be between 0 and 2pi and with a positive amplitude. 
% Options allow the phase to be between -pi and pi, with a positive amplitude, 
% or between 0 and pi with potentially a negavitve amplitude.
%
% syntax: [phase, amplitude] = Sinusoid_tidy(phase, amplitude, options)
%  options: 'PlusMinusPi'       -- force pi to be between -pi and +pi
%           'NegativeAmplitude' -- force phase to be 0 to pi. 
% phase and ampltude can be a single number or an array of numbers.
%
% Simon Hunt  December 2015. 


function [phases, amp] = Sinusoid_tidy(phases, amp, varargin)

    %% settings
    %default settings
    pm = 0;
    lim_low = 0;
    lim_high = 2*pi;
    pos_amp = 1;

    % Process the optional arguments overiding defaults
    iarg = 1 ;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'plusminuspi'
                pm = 1;
                lim_low = -pi;
                lim_high = pi;
                iarg = iarg + 1;
            case 'negativeamplitude'
                pos_amp = 0;
                iarg = iarg + 1;
            otherwise
                error(['Unknown option: ' varargin{iarg}]) ;
        end
    end
    
    %% input validation
    if size(phases) ~= size(amp)
        error('The number of phases and amplitudes do not match')
    end
    
    
    %% check amplitudes.
    % force all the ampltudes to be positive -- if required change back later.
    
    %find amplitudes less than 0
    neg = amp < 0;
    
    % correct amplitudes
    amp(neg) = -amp(neg);
    % correct corresponding phases
    phases(neg) = phases(neg)+pi;

    
    %% Limit phases
    % force phases to be between 0 and 2*pi.
    phases = phases - floor(phases./ (pi*2))* (2*pi);
    
    
    %% process input options: PlusMinusPi / negative amplitude    
    
    %find phases outside limits.
    large = phases >= pi;
    small = phases < 0;
        

    
    if pos_amp == 0 % if phase need to be 0<x<pi
        
        phases(large) = phases(large) -pi;
        amp(large)    = -amp(large);
        
        phases(small) = phases(small) +pi;
        amp(small)    = -amp(small);
        
    
    elseif pm == 1 %if phase to be +/-pi
        
        phases(large) = phases(large) -2*pi;
    
    end

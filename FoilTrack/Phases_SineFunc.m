function [ sine_array ] = Phases_SineFunc( time_array, period, amp1, varargin)
%Phases_SineFunc. Base sinusiod function for Phases Cancluaiton 
%
%   Syntax:
%     [values] = Phases_SineFunc(TimeStamps, period, A1, [P1], [A2], [P2])
%         values:       array of calculated sine values
%         TimeStamps:   array of time stamps
%
%         period:        period of sine function.
%         A1:            amplitude of the sine wave.
%         P1 (optional): phase offset of the sine wave (if absent assumed to be 0).
%         A2 (optional): amplitude of a sine^2 component (if absent assumed to be 0). 
%         P2 (optional): phase offset of the sine wave..
%  -  If there are three inputs a series of values of size(TimeStamps) is
%     returned assuming a phase of 0.
%  -  If four inputs the sine function has phase P1
%  -  If five inputs an additional sine^2 term is added to the function
%  with phase P1 and amplitude A2
%  -  If six inputs the sine^2 term has phase P2.
% 
%   See also: ImageAnalysis, PhaseCheck, Phases, PhasesSSD
   
% Simon Hunt. July 2017.
%  Version 1
          
%get size of time_stamps
[ht,wt] = size(time_array);
nt = numel(time_array);
if ~isvector(time_array)
    error('The time stamp array is not a vector')
end
time_array = time_array(:);

%get size of amplitude array
na = numel(amp1);

%parse inputs
if nargin >= 4
    phase1 = varargin{1};
    if numel(phase1) ~= numel(amp1)
        error('The number of phases does not match the number of amplitudes')
    end
else
    phase1 = zeros(1,na);
end
if nargin >= 5
    amp2 = varargin{2};
    if numel(amp2) ~= numel(amp1)
        error('The number of sine sq amplitudes does not match the number of amplitudes')
    end
else 
    amp2=zeros(1,na);
end
if isequal(amp2,zeros(1,na)) 
    phase2 = zeros(1,na);
elseif nargin >= 6
    phase2 = varargin{3};
    if numel(phase2) ~= numel(amp1)
        error('The number of sine sq phases does not match the number of amplitudes')
    end
else
    phase2=phase1;
end

%make sure arrays are all the right shape.
amp1 = amp1(:)';
phase1 = phase1(:)';
amp2 = amp2(:)';
phase2 = phase2(:)';

%calculate sine values.

%The sinusiod function in Phases.m was:
%   Y = amplitude(1) .* sin(X - phase) + amplitude(2) .* cos(2*(X - phase));
%   where cos(2t) was deemed to be equivalent to sin^2(t). It is like this
%   because the heating at X17B2 was current controlled and sinusoidal.
%   Temperature prop Power = I^2.R
%       I = A + a1.sin(tP) where P is the period.
%       I^2 => A + a'.sin(tP) + a''.sin^2(tP) prop T
%   In Phases.m the sin^2 term was replaced with 
%           => A + a'.sin(tP) + a''.cos(2tP)
%   Which is wrong becuase the cos is out of phase! It should have been
%   replaced with the proper sine identity:
%       sin^2(t) = 1/2(1-cos(2t)).
%
%whilst in PhasesSSD.m it was:
%   disp_min = amplitude .* cos(times./period*2*pi + phase);
%       there was never a need for the sin^2 term in the SSD data -- at
%       least it was never implemented. 
%
% There is no good reason for the differences between the equations or the stupid 
% errors in them (hence this script) -- at least the errors are now all the
% same.
%
% - I will stick to sine rather than cosine becuase most of the data sets are
% driven by sine rather than cos functions. 
% - Adding the phase offset rather than subtracting seems sensible at the time of writing.
% Especially because Kappa_solve has a line phases = -phases in it. 

% calculate the sine function:
% sine_array1 = repmat(amp1,nt,1).*sin(repmat(time_array./period*2*pi,1,na) + repmat(phase1,nt,1)) +...
%     repmat(amp2,nt,1).*1/2.*(1-cos(2*(repmat(time_array./period*2*pi,1,na) + repmat(phase2,nt,1))));

%this is a more computationally efficient copy of the equation above.
sine_array = amp1(ones(nt,1),:).*sin(time_array(:,ones(1,na))./period*2*pi +  phase1(ones(nt,1),:)) +...
    amp2(ones(nt,1),:).*1/2.*(1-cos(2*(time_array(:,ones(1,na))./period*2*pi + phase2(ones(nt,1),:))));

%check array has the same orientation as original time stamp array.
[hs,ws] = size(sine_array);
if ht == nt && hs~= nt
    sine_array = sine_array';
elseif wt == nt && ws~= nt
    sine_array = sine_array';
end









        
        
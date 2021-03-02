function [ out_amp, out_phase, varargout] = Sinusoid_Difference(amp1, phase1, amp2, phase2, varargin) 
% computes difference between two sine waves. 
% Returns single phase and amplitude for the function y = a sin(x+A) + b sin(x+B)
% assumes that the angles are in radians.
% after http://dspguru.com/sites/dspguru/files/Sum_of_Two_Sinusoids.pdf
% in \Simon\WORK folder (correct Nov 2014).
% syntax: 
%     [ out_amp, out_phase, varargout] = ...
%             Sinusoid_Difference(amp1, phase1, amp2, phase2, varargin) 
% Simon Hunt 2014

if nargin == 8
    [ out_amp, out_phase, out_amp_err, out_phas_err] = Sinusoid_Combine(amp1, phase1, -amp2, phase2, varargin{:}) ;
    varargout = {out_amp_err, out_phas_err};
else
    [ out_amp, out_phase] = Sinusoid_Combine(amp1, phase1, -amp2, phase2);
end

% x = 0:pi/2:2*pi;
% 
% In1= amp1.*sin(x+phase1);
% In2= amp2.*sin(x+phase2);
% out = out_amp .*sin(x+out_phase);
% 
% diffs = abs(out - (In1-In2));
% largest = max(diffs);
% if largest > out_amp/1e4
%     out_phase = out_phase+pi;
% end


if (0)
    x = 1:pi/50:2*pi;
    In1= amp1.*sin(x+phase1);
    In2= amp2.*sin(x+phase2);
    
    out = out_amp .*sin(x+out_phase);
    
    plot(x,In1,'g', x,In2,'g:', x,out,'b', x, In1-In2,'r:')
    
end

% out(15) - (In1(15)-In2(15));
% if abs(out(15) - (In1(15)-In2(15))) > 1e-4
%     keyboard
% end




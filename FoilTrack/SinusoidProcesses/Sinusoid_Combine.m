function [ out_amp, out_phase, varargout] = Sinusoid_Combine(amp1, phase1, amp2, phase2, varargin) 
% computes single phase and amplitude for the function y = a sin(x+A) + b sin(x+B)
% assumes that the angles are in radians.
% functions are taken from: http://dspguru.com/sites/dspguru/files/Sum_of_Two_Sinusoids.pdf
% the pdf is in \Simon\WORK folder (correct Nov 2014).
%
% version 2.0.2, January 2016
% Simon Hunt 2014-16


% preprocess inputs.
% need to make sure that the amplitudes are positive otherwise the answer
% is definitly wrong. 
[phase1, amp1] = Sinusoid_tidy(phase1, amp1);
[phase2, amp2] = Sinusoid_tidy(phase2, amp2);


%catch special case of sinusoides with same amplitude that are out of phase 
% should work within following functions but often fails due to
% numberical precision errors. Therefore have error catch.
if amp1 == -amp2 || (amp1 == amp2 && isclose(abs(phase1 -phase2), pi))
    out_amp = 0;
    out_phase = NaN;
    if nargin == 8
        varargout = {NaN(1) NaN(1)};
    end
    return
end

if exist([matlabroot '/toolbox/symbolic'], 'dir') == 7
    
    %set the equations using symbobolic logic toolbox.
    syms a1 p1 a2 p2
    amp = sqrt((a1 .*cos(p1) + a2.*cos(p2)).^2 + (a1 .*sin(p1) + a2.*sin(p2)).^2 );
    phase = atan( (a1.*sin(p1) + a2.*sin(p2)) ./ (a1.*cos(p1) + a2.*cos(p2)) );
    
    varlist = [a1 p1 a2 p2];
    vals = [amp1 phase1 amp2 phase2];
    out_amp = subs(amp,varlist,vals);
    out_phase = subs(phase,varlist,vals);
   
    if isa(out_amp,'sym') == 1
        out_amp   = eval(out_amp);
        out_phase = eval(out_phase);
    end
        
else
    
    % compute combined sine / cos function
    out_amp = sqrt((amp1 .*cos(phase1) + amp2.*cos(phase2)).^2 + (amp1 .*sin(phase1) + amp2.*sin(phase2)).^2 );
    out_phase = atan( (amp1.*sin(phase1) + amp2.*sin(phase2)) ./ (amp1.*cos(phase1) + amp2.*cos(phase2)) );
    
end


%check combined function has correct phase (does not suffer from phase wrapping)
x = 0:pi/10:2*pi;

In1= amp1.*sin(x+phase1);
In2= amp2.*sin(x+phase2);
out = out_amp .*sin(x+out_phase);

diffs = abs(out - (In1+In2));
largest = max(diffs);
if largest > out_amp/1e4
    out_phase = out_phase+pi;
    corr_phase = 1;
else
    corr_phase = 0;
end

%plot input and output sinusoids. Used for debugging.
bebug = 0;
if bebug == 1
    out = out_amp .*sin(x+out_phase);
    plot(x,In1,'g', x,In2,'g:', x,out+0.01,'b', x, In1+In2+0.02,'r')
    legend('In 1', 'In 2', 'combined', 'sum')
%     keyboard
end



%error propagation
if ~isempty(varargin)
    if exist([matlabroot '/toolbox/symbolic'], 'dir') == 7
        %calculate errors from symbolic logic.
        OutAmpErr   = PropError(amp,varlist,vals,cell2mat(varargin));
        OutPhaseErr = PropError(phase,varlist,vals,cell2mat(varargin));
    else
        % done numerically by monte carlo.
        
        amp1err = varargin{1};
        phase1err = varargin{2};
        amp2err = varargin{3};
        phase2err = varargin{4};
        
        number = 1E7; % number of values to use in computing the errors
        
        %make arrays of values
        amp1_array = randn(number,1)*amp1err + amp1;
        phase1_array = randn(number,1)*phase1err + phase1;
        amp2_array = randn(number,1)*amp2err + amp2;
        phase2_array = randn(number,1)*phase2err + phase2;
        
        %compute array of phases and amplitudes
        OutAmpArray = sqrt((amp1_array .*cos(phase1_array) + amp2_array.*cos(phase2_array)).^2 + (amp1_array .*sin(phase1_array) + amp2_array.*sin(phase2_array)).^2 );
        OutPhaseArray = atan( (amp1_array.*sin(phase1_array) + amp2_array.*sin(phase2_array)) ./ (amp1_array.*cos(phase1_array) + amp2_array.*cos(phase2_array)) );
        %catch incorrect error when wrapping around 0 degrees. 
        %not suite sure why it works
        OutPhaseArray(OutPhaseArray<=pi/2) = OutPhaseArray(OutPhaseArray<=pi/2)+pi;
        %end
        if corr_phase == 1
            OutPhaseArray(OutPhaseArray<0) = OutPhaseArray(OutPhaseArray<0) +pi;
        end
        
        %standard deviation of values. 
        OutAmpErr = std(OutAmpArray);
        OutPhaseErr = std(OutPhaseArray);
        
        %     figure(3), hist(OutAmpArray,20)
        %     figure(4), hist(OutPhaseArray,20)        
        %     figure(6), plot(OutAmpArray, OutPhaseArray, 'r.')
    
    end
    varargout = {OutAmpErr, OutPhaseErr};
end


end

function out = isclose(value, ideal_value)
    out = (abs(value-ideal_value) < 1e4*eps(min(abs(value),abs(ideal_value))));
end


function error = PropError(f,varlist,vals,errs)
%SIGMA = PROPERROR(F,VARLIST,VALS,ERRS)
%
%Finds the propagated uncertainty in a function f with estimated variables
%"vals" with corresponding uncertainties "errs".
%
%varlist is a row vector of variable names. Enter in the estimated values
%in "vals" and their associated errors in "errs" at positions corresponding 
%to the order you typed in the variables in varlist.
%
%Example using period of a simple harmonic pendulum:
%
%For this example, lets say the pendulum length is 10m with an uncertainty
%of 1mm, and no error in g.
%syms L g
%T = 2*pi*sqrt(L/g)
%type the function T = 2*pi*sqrt(L/g)
%
%PropError(T,[L g],[10 9.81],[0.001 0])
%ans =
%
%    [       6.3437]    '+/-'    [3.1719e-004]
%    'Percent Error'    '+/-'    [     0.0050]
%
%(c) Brad Ridder 2007. Feel free to use this under the BSD guidelines. If
%you wish to add to this program, just leave my name and add yours to it.
n = numel(varlist);
sig = vpa(ones(1,n));
for i = 1:n
    sig(i) = diff(f,varlist(i),1);
end
error1 =sqrt((sum((subs(sig,varlist,vals).^2).*(errs.^2))));
error = double(error1);
% sigma = [{subs(f,varlist,vals)} {'+/-'} {error};
%          {'Percent Error'} {'+/-'} {abs(100*(error)/subs(f,varlist,vals))}];
end



%%versions
% -v2.0.2  13 Jan 2016
%   - made number of arguments out match the number of arguments incase of
%   amplitudes being equal and oposit. 
% -v2.0.1  13 Jan 2016
%   - minor changes to make script work with R2104a
% -v2  Jan 2016
%   - rewote the maths using the symbolic logic toolbox. This is
%   approximately an order of magnitude faster than the monte carlo verion.
% -v1  
%   - combine two sines using equaitons from pdf listed in header.
%   - errors calculated using monte carlo method
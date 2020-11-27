function Sinusoid_tidy_TestCases(varargin)
%SinusoidProcesses_TestCases -- test outputs of all SinusoidProcess funtions 
% checks that the outputs of Sinusoid_Tidy are within numberical error of
% the expected results.
%       Simon Hunt 31st December 2015

test_num = 0;

%% Sinusoid_Tidy

%1. single positive phase <pi and positive amplitude
test_num = test_num+1;
phase = 1;
amp = .5;
fprintf('Test %d: Sinusiod_tidy(0<phase<pi, amplitude>0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi && min(phase_out(:)) >= 0 ...
        && phase_out == phase && amp_out == amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end

%2. single positive phase pi<x<2pi and positive amplitude
test_num = test_num+1;
phase = 4;
amp = .5;
fprintf('Test %d: Sinusiod_tidy(0<phase<pi, amplitude>0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi && min(phase_out(:)) >= 0 ...
        && phase_out == phase && amp_out == amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%3. single phase (0<x< pi) with negtive amplitude
test_num = test_num+1;
phase = 2;
amp = -.5;
fprintf('Test %d: Sinusiod_tidy(0<phase<pi, amplitude>0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi && min(phase_out(:)) >= 0 ...
        && phase_out == phase+pi && amp_out == -amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%4. single phase (pi<x<2pi) with negative ampltude
test_num = test_num+1;
phase = 6;
amp = -.5;
fprintf('Test %d: Sinusiod_tidy(0<phase<pi, amplitude>0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi && min(phase_out(:)) >= 0 ...
        && phase_out == phase+pi-2*pi && amp_out == -amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%5. single phase x>>2pi and positive amplitude
test_num = test_num+1;
phase = .1+12*pi;
amp = .05;
fprintf('Test %d: Sinusiod_tidy(phase>>2pi, amplitude>0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi && min(phase_out(:)) >= 0 ...
        && phase_out == phase - 12*pi && amp_out == amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%6. single phase x>>2pi and negative amplitude
test_num = test_num+1;
phase = .1+13*pi;
amp = -.05;
fprintf('Test %d: Sinusiod_tidy(phase>>2pi, amplitude<0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi && min(phase_out(:)) >= 0 ...
        && phase_out == phase+pi-7*2*pi && amp_out == -amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%7. single phase x<<0 and positive amplitude
test_num = test_num+1;
phase = .1-12*pi;
amp = .05;
fprintf('Test %d: Sinusiod_tidy(phase<<0, amplitude>0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if phase_out == phase+12*pi && amp_out == amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%8. single phase x<<0 and negative amplitude
test_num = test_num+1;
phase = .1-12*pi;
amp = -.05;
fprintf('Test %d: Sinusiod_tidy(phase<<0, amplitude<0)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi && min(phase_out(:)) >= 0 ...
        && isclose(phase_out, .1+pi) && amp_out == -amp
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%9. Array of phase_out and ampltudes.
test_num = test_num+1;
phase = [.1-12*pi 1 2 12; -6 0 -1 .5];
amp = [-.05 1 2 -.1; 1 4 0 0.002];
phase_check = [.1+pi 1 2 12+pi-4*pi; -6+2*pi 0 -1+2*pi .5];
amp_check   = [.05   1 2 .1;         1      4  0      0.002];
fprintf('Test %d: Sinusiod_tidy(phase_array, amplitude_array)\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp);
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < 2*pi & min(phase_out(:)) >= 0 ...
        & isclose(phase_out, phase_check) & amp_out == amp_check
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Tidy does not return same numbers when run under normal contidions')
end


%10. plusminus pi options
test_num = test_num+1;
phase       = [.1-12*pi 1 2 12;         -6      0 -1 .5];
phase_check = [.1-pi    1 2 12+pi-4*pi; -6+2*pi 0 -1 .5];
amp         = [-.05     1 2 -.1;        1       4  0 0.002];
amp_check   = [.05      1 2 .1;         1       4  0 0.002];
fprintf('Test %d: Sinusiod_tidy(phase_array, amplitude_array, ''PlusMinusPi'')\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp, 'PlusMinusPi');
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < pi & min(phase_out(:)) >= -pi ...
        & isclose(phase_out, phase_check) & amp_out == amp_check
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    phase, phase_out, phase_check, amp, amp_out
    error('Sinusoid_Tidy does not return same numbers when run with the option ''PlusMinusPi''')
end


%10. negative amplitude options
test_num = test_num+1;
phase       = [.1-12*pi 1+pi; 2 12;         -6      0; -1 .5];
phase_check = [.1       1; 2 12-3*pi; -6+2*pi 0; -1+pi .5];
amp         = [-.05     1; 2 -.1;        1       -4;  0 0.002];
amp_check   = [-.05     -1; 2 .1;         1       -4;  0 0.002];
fprintf('Test %d: Sinusiod_tidy(phase_array, amplitude_array, ''NegativeAmplitude'')\n', test_num)
fprintf('   Input:    phase = %4.2f and amplitude = %4.2f\n', phase, amp)
[phase_out amp_out] = Sinusoid_tidy(phase, amp, 'NegativeAmplitude');
fprintf('   Returned: phase = %4.2f and amplitude = %4.2f\n', phase_out, amp_out)
if max(phase_out(:)) < pi & min(phase_out(:)) >= 0 ...
        & isclose(phase_out, phase_check) & amp_out == amp_check
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    phase, phase_out, phase_check, amp, amp_out
    error('Sinusoid_Tidy does not return same numbers when run with the option ''NegativeAmplitude''')
end




fprintf('\n\n   Sinusoid_tidy: All tests cases passed.\n\n')




end

function out = isclose(value, ideal_value)
    out = (abs(value-ideal_value) < 1e4*eps(min(abs(value),abs(ideal_value))));
end



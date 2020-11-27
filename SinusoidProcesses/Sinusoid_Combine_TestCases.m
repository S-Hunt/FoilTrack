function Sinusoid_Combine_TestCases(varargin)
%Sinusoid_Combine_TestCases -- test output of Sinusoid_combineSinewaves  
% checks that the outputs of Sinusoid_Tidy are within numberical error of
% the expected results.
%       Simon Hunt 31st December 2015

if ~isempty(varargin)
    test_function = str2func(varargin{1});
else
    test_function = @Sinusoid_Combine;
end

test_num = 0;
    
%1. sine wave + sine with 0 amplitude
test_num = test_num+1;
phase1 = 1.1234;
amp1 = 0.024;
phase2 = 0;
amp2 = 0;
fprintf('Test %d: Sine wave with sine wave of no amplitude\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(phase_out, phase1) & isclose(amp_out, amp1)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause


%2. sine wave + sine with 0 amplitude (same as test 1 but other way round)
test_num = test_num+1;
phase2 = rand(1);
amp2 = rand(1);
phase1 = 0;
amp1 = 0;
fprintf('Test %d: Sine wave with sine wave of no amplitude, other possible variation\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(phase_out, phase2) & isclose(amp_out, amp2)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end

%pause

%3. sine wave (phase>>2pi) + sine with 0 amplitude
test_num = test_num+1;
phase2 = 12*pi;
amp2 = rand(1);
phase1 = 0;
amp1 = 0;
fprintf('Test %d: Sine wave (phase>>2pi) with sine wave of no amplitude\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(phase_out, 0) & isclose(amp_out, amp2)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause


%4. sine wave (amp<0) + sine with 0 amplitude
test_num = test_num+1;
phase2 = rand(1);
amp2 = -rand(1);
phase1 = 0;
amp1 = 0;
fprintf('Test %d: Sine wave of no amplitude with sine wave with -ve amplitude\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(phase_out, phase2+pi) & isclose(amp_out, -amp2)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    fprintf('\n TEST FAILED. \n\n')
    fprintf('   Returned values:    phase1 = %4.2f,  amp1 = %4.2f\n', phase_out, amp_out)
    fprintf('   Correct values:     phase1 = %4.2f,  amp1 = %4.2f\n', phase2+pi, -amp2)
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause


%5. sine wave + sine with same amplitude 
test_num = test_num+1;
phase1 = rand(1);
amp1 = rand(1);
phase2 = phase1;
amp2 = amp1;
fprintf('Test %d: Sine wave combined with itself\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(phase_out, phase1) & isclose(amp_out, amp1*2)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause



%6. sine wave - sine with same amplitude 
test_num = test_num+1;
phase1 = rand(1);
amp1 = rand(1);
phase2 = phase1;
amp2 = -amp1;
fprintf('Test %d: Sine wave combined with itself (with negative amplitude)\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(amp_out, 0)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause

%7. sine wave + sine with same amplitude out of phase by pi
test_num = test_num+1;
phase1 = rand(1);
amp1 = rand(1);
phase2 = phase1+pi;
amp2 = amp1;
fprintf('Test %d: Sine wave combined with wave of sample amplitude out of phase by pi\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(amp_out, 0)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause


%8. sine wave + sine with same amplitude out of phase by pi
test_num = test_num+1;
phase1 = rand(1);
amp1 = rand(1);
phase2 = phase1+pi/2;
amp2 = amp1;
fprintf('Test %d: Sine wave combined with wave of sample amplitude out of phase by pi/2\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(amp_out, amp1*sqrt(2)) & isclose(phase_out, phase1+pi/4)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause



%9. sine wave + sine with half -ve amplitude
test_num = test_num+1;
phase1 = rand(1);
amp1 = rand(1);
phase2 = phase1;
amp2 = -amp1/2;
fprintf('Test %d: Sine wave combined with wave of negative half amplitude in phase\n', test_num)
[amp_out phase_out] = test_scenario(test_function, amp1, phase1, amp2, phase2);
if isclose(amp_out, amp1/2) & isclose(phase_out, phase1)
    fprintf(' Passed test %d: input and output values correspond\n\n', test_num)
else
    error('Sinusoid_Combine does not return the expected values when run')
end
%pause



fprintf('\n\n   Sinusoid_Combine: All tests cases passed.\n\n')




end


function [a_out p_out] = test_scenario(f, a1, p1, a2, p2)

fprintf('   %s(amp1, phase1, amp2, phase2)\n', func2str(f))
fprintf('   Input:    phase1 = %4.2f, amp1 = %4.2f\n', p1, a1)
fprintf('             phase2 = %4.2f, amp2 = %4.2f\n', p2, a2)
[a_out p_out] = f(a1, p1, a2, p2);
fprintf('   Returned: phase = %4.2f,  amplitude = %4.2f\n', p_out, a_out)

end


function out = isclose(value, ideal_value)
    out = (abs(value-ideal_value) < 1e4*eps(min(abs(value),abs(ideal_value))));
end



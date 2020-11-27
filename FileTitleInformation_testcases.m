%% FileTitleInformation_testcases
%
% This script runs a series of testcases to get information from strings used in data processing. 
%   It must always work with all previous cases.
% 
% See Also: FileTitleInformation

%   Simon Hunt 2018
%       version 1




function [out] = FileTitleInformation_testcases



%% Case 1: Typical NSLS/APS file name 
pass = 0;
text_str = 'Garn43_20pt4tons_300C_10s_004.nc';
disp('Case 1: Typical NSLS/APS file name');


t = FileTitleInformation(text_str);

if t{1} == 'Garn43' &...
    t{2} == 20.4 &...
    t{3} == 300 &...
    t{4} == 10 &...
    t{5} == 4

pass = 1; 
end

if pass ~= 1
    keyboard
end


%% Case 2: difficult NSLS/APS file name 
pass = 0;
text_str = 'CaI_20_20pt4tons_1000C_10s_324_sine_fits.txt';
disp('Case 1: difficult NSLS/APS file name');


t = FileTitleInformation(text_str);

if t{1} == 'CaI_20' &...
    t{2} == 20.4 &...
    t{3} == 1000 &...
    t{4} == 10 &...
    t{5} == 324

pass = 1; 
end

if pass ~= 1
    keyboard
end



%% Case 3: Thermal expansion of Ice name 
pass = 0;
text_str = 'D2O_IceIh_20K.lst';
disp('Case 1: difficult thermal expansion name');


t = FileTitleInformation(text_str);

if strcmpi(t{1}, 'D2O_IceIh') &...
    isnan(t{2}) &...
    t{3} == 20-273 &...
    isnan(t{4}) &...
    isnan(t{5})

pass = 1; 
end

if pass ~= 1
    keyboard
end




%% Case 4: dirctory names
pass = 0;
text_str = '100C';
disp('Case 4: directory name just temperture (100C)');


t = FileTitleInformation(text_str);

if strcmpi(t{1}, 'Not given') &...
    isnan(t{2}) &...
    t{3} == 100 &...
    isnan(t{4}) &...
    isnan(t{5})

pass = 1; 
end

if pass ~= 1
    keyboard
end


%% Case 5: dirctory names
pass = 0;
text_str = '100 Kelvin';
disp('Case 5: directory name just temperture (100 kelvin)');


t = FileTitleInformation(text_str);

if strcmpi(t{1}, 'Not given') &...
    isnan(t{2}) &...
    t{3} == 100-273 &...
    isnan(t{4}) &...
    isnan(t{5})

pass = 1; 
end

if pass ~= 1
    keyboard
end

%% Case 6: dirctory names
pass = 0;
text_str = '100_Kelvin';
disp('Case 6: directory name just temperture (100_kelvin)');


t = FileTitleInformation(text_str);

if strcmpi(t{1}, 'Not given') &...
    isnan(t{2}) &...
    t{3} == 100-273 &...
    isnan(t{4}) &...
    isnan(t{5})

pass = 1; 
end

if pass ~= 1
    keyboard
end


disp('All tests passed')



% file_name_cut = strrep(file_name_cut,'_kelvin','k');
% file_name_cut = strrep(file_name_cut,'kelvin','k');


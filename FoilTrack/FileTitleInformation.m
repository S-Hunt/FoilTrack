%% FileTitleInformation
%
% This script reads the file names of experiments from the NSLS or other locations and returns
% the parameters listed in the file name.
%
% [experiment_name load temperature|power period (file number in series)...
%       (other information) (file name ending)] = FileTitleInformation('filename');
% where 'filename' is of the type Exp_003_40T500C1s_right.nc or similar
%
% This is done by finding the numbers preceeding the 'T'/'bars', 'C'/'W' and
% 's'/'Hz' in the file name.
% 
% See Also: ImageAnalysis

%   Simon Hunt 2009 - 2018
%       version 3.5.2

% Known bugs: only deals with files without extra stuff in titles to experimental information (i.e. *.nc files) or txt files.


function [out prop_array] = FileTitleInformation(file_name, varargin)


% Process the input arguments
iarg = 1 ;
while iarg <= (length(varargin))
    switch varargin{iarg}
        case 'warn off'
            s.a = warning('off', 'FileTile:noTemp');
            s.b = warning('off', 'FileTile:noLoad');
            s.c = warning('off', 'FileTile:noPeriod');
            iarg = iarg+1;
        otherwise
%         keyboard
        error(['Unknown option: ' varargin{iarg}]) ;
    end

end

%remove leading or traling white space
file_name = strtrim(file_name);

%replace all spaces with '_'
file_name = strrep(file_name, ' ', '_');

% split constituent parts of the file name into parts.
[~, file_name_str_Upper, type] = fileparts(file_name);
file_name_str = lower(file_name_str_Upper);
type = strtrim(type(2:end));
f_type = [];

tpye_long = strfind(type,',');
if ~isempty(tpye_long)
    type = type(1:tpye_long-1);
end

switch lower(type)
    case {'txt' 'mat' 'tif'}
        %remove extra titles added by Phases and ImageAnalysis
        titl = min([strfind(file_name_str,'sine_fits') ...
                     strfind(file_name_str,'position_change')...
                     strfind(file_name_str,'sine_sq_fits')...
                     strfind(file_name_str,'boxes')...
                     strfind(file_name_str,'ssd')]);
        
        if titl ~= 1
            file_name_cut = file_name_str(1:titl-1);
            while strcmpi(file_name_cut(end),'_') == 1
                file_name_cut = file_name_cut(1:end-1);
            end
        elseif file_name_str(1) == 's' %this catch was for some expeirment -- unknown now
            file_name_cut = file_name_str(10:end);
        elseif file_name_str(1) == 'p' %this catch was for some expeirment -- unknown now
            file_name_cut = file_name_str(16:end);
        else
            file_name_cut = file_name_str;
%             error('I have tried to remove extra text from file title and failed')
        end

        if file_name_cut(1) == '_'
            file_name_cut = file_name_cut(2:end);
        end
        
        %make string for file type output.
        f_type = file_name_str_Upper(titl:end);
        
    case 'nc'
        file_name_cut = file_name_str;
        %this is otherwise left blank delibarately.
        
    case 'lst'
        file_name_cut = file_name_str;
        %this is otherwise left blank delibarately.

    otherwise
        file_name_cut = file_name_str;
%         warning('The file type is unrecognised')
%         out = [{file_name_str},{NaN},{NaN},{NaN},{NaN}];
%         return
end


% replace all pt in names with point (used for non integer values)
%assume that it is only a value if lead and trailed by a number
pts = strfind(file_name_cut,'pt');
if ~isempty(pts)
    for x = 1:length(pts)
        %check if character before and after the 'pt' are numbers.
        stend(x) = isempty(str2num(file_name_cut(pts(x)-1))) + isempty(str2num(file_name_cut(pts(x)+2)));
    end
    remove = (stend ~= 0);
    pts(remove) = [];
end
%replace the pt surrouneded by numbers with '.'
for x = length(pts):-1:1
    file_name_cut(pts(x)+1) = [];
    file_name_cut(pts(x)) = '.';
end
% file_name_cut = strrep(file_name_cut,'pt','.');

%replace all the kelvins with K.
file_name_cut = strrep(file_name_cut,'_kelvin','k');
file_name_cut = strrep(file_name_cut,'kelvin','k');



%categorise characters in names into numbers and letters.
% all '.' are characterised as numbers to allow for recongitiion of
% non-integer values.
for j = 1:length(file_name_cut)
    which(j) =  isdigit(file_name_cut(j));
    if file_name_cut(j) == '.'
        which(j) = 1;
    end
end

numbers = find(which == 1);
letters = find(which == 0);


%% Temperature of the files

%find location of C and W for regonition of where the temperature is
degC_pos = find((file_name_cut == 'c'),1,'last');
K_pos = find(file_name_cut == 'k');
W_pos = find(file_name_cut == 'w');
RT_pos =  strfind(file_name_cut, 'rt');

%if isempty(degC_pos) == 0 && (degC_pos >  find(file_name_cut=='_',1,'first') ) && isempty(str2num(file_name_cut(degC_pos-1)))==0 %...%if there is a celcius temperature
if ~isempty(degC_pos)...
        && degC_pos(1) ~= 1 ... the C is not the first character in the name - needed to chatch pos-1 searchers
        && isempty(str2num(file_name_cut(degC_pos-1)))==0 ...
        && strcmpi(file_name_cut(degC_pos-1),'i') == 0 %% i is a number so have to check that the value in front is not an i. (this occurs for example in 'ice')
    %&& isnumeric(file_name_cut(degC_pos-1))==1 %if numbers are in front of the c
    start_temp = letters(find(letters < degC_pos,1,'last'));
    if isempty(start_temp)
        start_temp = 0;
    end
    if start_temp ~= degC_pos - 1
        temp = str2num(file_name_cut(start_temp+1 : degC_pos-1));
        prop_array.temperature = str2num(file_name_cut(start_temp+1 : degC_pos-1));
        prop_array.temperature_units = 'degC';
        typeT(1) = 1;
    else
        start_temp = NaN;
        temp = NaN;
        typeT = [];
    end
elseif isempty(K_pos) == 0  && isempty(str2num(file_name_cut(K_pos-1)))==0 %if the temperature is in kelvin
    start_temp = letters(find(letters < K_pos,1,'last'));
    if isempty(start_temp)
        start_temp = 0;
    end
    temp = str2num(file_name_cut(start_temp+1 : K_pos-1))-273; %return temp in degC
    prop_array.temperature = str2num(file_name_cut(start_temp+1 : K_pos-1));
    prop_array.temperature_units = 'K'; %store temperature in array in Kelvin.
    typeT(3) = 1;
elseif ~isempty(W_pos)... %if the temperature is a power
        && W_pos(1) ~= 1 ... the C is not the first character in the name - needed to chatch pos-1 searchers
        && isempty(str2num(file_name_cut(W_pos-1)))==0 ...
        && strcmpi(file_name_cut(W_pos-1),'i') == 0 %% i is a number so have to check that the value in front is not an i. (this occurs for example in 'ice')
    start_temp = letters(find(letters < W_pos,1,'last'));
    temp = str2num(file_name_cut(start_temp+1 : W_pos-1));
    prop_array.temperature = str2num(file_name_cut(start_temp+1 : W_pos-1));
    prop_array.temperature_units = 'W';    
    typeT(2) = 1;
elseif isempty(RT_pos) == 0 %if the temperature is written as RT
    start_temp = RT_pos;
    temp = 22;
    prop_array.temperature = 22;
    prop_array.temperature_units = 'degC';    
    typeT(2) = 1;
else
    warn_msg = ['No recognised temperatures in the file name: ',file_name_str_Upper, '.', type];
    warning('FileTile:noTemp',warn_msg);
    start_temp = NaN;
    temp = NaN;
    typeT = [];
end
if sum(typeT) == 2 %displays a warning if temp and power are mixed in the experiment files.
    warning('There is a combination of power and/or temperature units in this experiment.');
    typeT(3) = 1;
end


%% Load/Force information in title

%find the string denoting the load
load_pos_GPa = strfind(file_name_cut,'gpa'); %find pressure in GPa
load_pos_tonne = strfind(file_name_cut, 'tonne'); %find load in ton(ne)s
    if isempty(load_pos_tonne)
        load_pos_ton = strfind(file_name_cut, 'ton');
    else
        load_pos_ton = [];
    end
load_pos_bar = strfind(file_name_cut,'bar'); %find load in bars
    if isempty(load_pos_ton) && isempty(load_pos_bar) && isempty(load_pos_tonne) && isempty(load_pos_GPa)
        load_pos_t = strfind(file_name_cut, 't');  
        if ~isempty(load_pos_t) && load_pos_t(1) == 1
            load_pos_t(1) = [];
        end
        for x = 1:length(load_pos_t)
            if ~isdigit(file_name_str(load_pos_t(x)-1)) && load_pos_t(x)~=1
                load_pos_t(x) = NaN;
            end
        end
        load_pos_t = load_pos_t(isfinite(load_pos_t));
    else 
        load_pos_t = [];
    end
load_pos = [load_pos_ton, load_pos_bar, load_pos_tonne, load_pos_t, load_pos_GPa];

%read information 
if ~isempty(load_pos)
    start_load = letters(find(letters < load_pos,1,'last'));
    if isempty(start_load)
        start_load = 0;
    end
    load = str2num(file_name_cut(start_load+1 : load_pos-1));
    prop_array.load = str2num(file_name_cut(start_load+1 : load_pos-1));
else
    warn_msg = ['No recognised load in the file name: ',file_name];
    warning('FileTile:noLoad',warn_msg);
    start_load = NaN;
    load = NaN;
end

%get units
if ~isempty(load_pos_tonne)
    prop_array.load_units = 'tonne';
elseif ~isempty(load_pos_ton)
    prop_array.load_units = 'ton';
elseif ~isempty(load_pos_bar)
    prop_array.load_units = 'bars';
elseif ~isempty(load_pos_GPa)
    prop_array.load_units = 'GPa';
end   

%% Period/Frequency information in title.

%find position of driving period (if it is present)
period_pos = strfind(file_name_cut, 's'); %find period in seconds
if ~isempty(period_pos)
    if isempty(period_pos)
        period_pos = strfind(file_name_cut, 'set_'); %this string appears in the CaI_* set of experiments.
    end
%     if period_pos(1) == 1 %If the experiment name starts with an 's'
%         period_pos(1) = [];
%     end

    %Ignore (remove) 's' from the list of posibilities if it comes before the first _
    und = strfind(file_name_cut, '_');
    if ~isempty(und)
        period_pos(period_pos < und(1)) = [];
    end
    if length(period_pos) >= 1 %if there is more than 1 s check the characters in front of it to see which is the period
        leading = file_name_cut(period_pos-1); %list of characters preceeding the s
        for x = 1 : length(period_pos)
            choose(x) = (isdigit(leading(x))...    %is the charater in front a number
                          | leading(x) == 'f'...   %is the character in front an f (for half)
                          | (leading(x) == 'r' &...%is the character in front an r (for quarter) but not the whole string is not 'bars'
                              strcmpi(file_name_cut(period_pos(x)-2:period_pos(x)-0),'ars')==0 ));
        end
        period_pos = period_pos(choose == 1);
    end
end
    
%get the period
if ~isempty(period_pos)
    start_period = letters(find(letters < period_pos,1,'last'));
    if period_pos - start_period == 1 %confirm that the periods are not half or quarter as strings. 
        if strcmpi(file_name_cut(period_pos-4:period_pos-1), 'half') == 1
            period = 0.5;
        elseif strcmpi(file_name_cut(period_pos-6:period_pos-1), 'quarter') == 1
            period = 0.25;
        else
            period = NaN;
        end
    else %the period should be a number...
            period_str = file_name_cut(start_period+1 : period_pos-1);
            period = str2num(period_str);
            if isempty(period) == 1
               error 'The period has not been found';
            end
    end
    prop_array.period = period;
    prop_array.period_units = 's';    
else %we think there is a period but cannot find it.
    warn_msg = ['No recognised periods in the file name: ',file_name];
    warning('FileTile:noPeriod',warn_msg);
    start_period = NaN;
    period = NaN;
end


%% experiment name
title_end = min([start_period start_load start_temp]);
if isnan(title_end)
    title_end = length(file_name_cut)+1;
end
expt_name = file_name_str_Upper(1:title_end-1);
if isempty(expt_name)
   expt_name = 'Not given'; 
end
prop_array.expt_name = expt_name;


%% run number in the title.
underscore = find(file_name_cut == '_');
if ~isnan(min([start_period start_load start_temp])) && ~isempty(underscore)
    RunNum = file_name_cut(underscore(end)+1:end);
    RunNum = str2double(RunNum);
    prop_array.run_number = RunNum;
elseif isnan([start_period start_load start_temp]) & ~isempty(underscore) %run name and number.
    RunNum = file_name_cut(underscore(end)+1:end);
    RunNum = str2double(RunNum);
    prop_array.run_number = RunNum;
else
    
    RunNum = NaN;
    prop_array.run_number = NaN;
end




%% extras
% this section is not as well defined and will evolve.

%find the position after which the end has to be.
conditions_end = underscore(find(underscore > max([start_period start_load start_temp])));

%get the extra text if it exists.
if contains(file_name_cut, 'right')
    extra = 'right';
elseif contains(file_name_cut, 'left')
    extra = 'left';
elseif contains(file_name_cut, 'middle')
    extra = 'middle';
elseif ~isempty(conditions_end) && isnan(RunNum) && conditions_end(1)==underscore(end) %if the RunNum is not a number read the end string in the file name
    extra = file_name_cut(underscore(end)+1:end);
elseif ~isempty(conditions_end) && conditions_end(1)==underscore(end-1) %if RunNum is a number read the second last string in the file name
    extra = file_name_cut(underscore(end-1)+1:underscore(end)-1);
elseif contains(file_name_cut, 'ref')
    s = strfind(file_name_cut, 'ref');
    p = [conditions_end, length(file_name_cut)];
    e = find(p > s, 1, 'first');
    extra = file_name_cut(s:p(e));
else
    extra = [];
end
if ~isempty(extra)
    prop_array.extra = extra;
end

%% constructs output array 
out = [expt_name,{load},{temp},{period},{RunNum},{extra},{f_type}];

% prop_array
end



function [out, org] = isdigit(in)

out = ~isempty(regexprep(in,'[^.0-9]',''));

end



%% version changes
% V 3.5.2   - 17th December 2018
%   - bug catch to allow experiments with names starting or containing 'W'.
% V 3.5.1   - 26th October July 2018
%   - bug catch to allow experiments with names starting in 'C'.
% V 3.5     - 3rd July 2018
%   - added lines to return ref number if presenet in the file name.
% V 3.4.5   - 11th May 2018
%   - replace all spaces in filename with underscores.
% V 3.4.4   - 7th May 2018
%   - minor changes to temperature finding so it works with 'ice'
% V 3.4.3   - 20th December 2017
%   - Changed temperature finding for degC condition to exclude 'c' from
%   first part of file name -- so that it works with Ice data from ISIS. I
%   replaced the logic condition '(degC_pos ~= 1)' with '(degC_pos >
%   find(file_name_cut=='_',1,'first') )'.
% V 3.4.2   - 17th August 2017
%   - Edited period getting lines so that the experiment name can validly
%   have an 's'.
% V 3.4.1   - 21st May 2017
%   - Added file type to the end of the output array. e.g. SSD, positon_change, etc.
% V 3.4   - 13th Feb 2017
%   - Added ability to read frequency in hertz. 
%   - Added ability to skip experiment name if it is not in the title.
% V 3.3.3 - 16th Jan 2017
%   - Added 'GPa' as a recognised ending to the load (now pressure as well).
%   - Replaced 'isdigit' function with a regularexpression statement to make it faster
% V 3.3.2 - 16th Jan 2017
%   - Bug fixes to do with 'pt' strings in file names.
% v 3.3.1 4 oct 2016
%     - added option to recognise RT at 25C.
% v 3.3.0 27th July 2016
%     - added option to supress the warnings generated when the files do
%     not contain all the expected title information.
% v 3.2.0 July 2016
%     - made the options for reading the extra information in the file name
%     generic.
% v 3.1.0 June 2016
%     - added kelvin options for temperature. 
% v 3.0.2 April 2016
%     - bug fixes for when names do not contain a full set of information.
% v 3.0.1 April 2016
%    - Adjusted the code so that it returns proper results from file names with no information in them
% v3 - version numbering started
%    - change script so that recognises pt and returns non integer values.
% WriteSineFit%
%   Writes the output of the sinusiods fitted to the displacements calculated by the 
%   ImageAnalysis script.
%
%   Syntax: WriteSineFit(OutFileName, variables_to_write, header_vairables_to_write)
%     variables_to_write is a structure which the fields are the variables to be saved
%       and its units. 
%     The variables saved are: 'period', 'amplitude', 'phase', 'horizontal_position_mm'
%     Optionally can take variables called: 'SSD', 'ChiSq', 'num_used_dat', 'use_in_phase_fit', 
%
%    See also: READSINEFIT, Phases, PhasesSSD

% Simon Hunt November 2014
% $ version 2.1 $ June 2016 $

function WriteSineFit(OutFileName, variables_to_write, header_vairables_to_write, varargin)


%% read values from caller function for header and create it. 

%read values
% file name data was read from.
filename = evalin('caller', 'file_name'); 

%get date of creation for data file.
input_file_details = dir(strcat(filename,'*'));
input_file_date = input_file_details.date;

%find name of function called to do analysis
RunningStack = dbstack('-completenames');
analysis_file = dir(RunningStack(2).file);

%reads version number and creation date from file doing analysis.
fid = fopen(analysis_file.name);
content = fscanf(fid,'%c', 1500);
dollars = find(content == '$');
Rev =  strtrim(content(dollars(1)+1:dollars(2)-1));
fclose(fid);

%add DEV to the output file name if the caller is a development version of
%the script (i.e. has dev in the file name or the version number). 
if (~isempty(strfind(lower(analysis_file.name), 'dev')) || ~isempty(strfind(lower(Rev), 'dev')) )...
    && isempty(strfind(lower(OutFileName), 'dev'))
    
    [location,root,ending] = fileparts(OutFileName);
    OutFileName = [location, root, '_dev', ending];
    
end


%% Output file header
options_list{1} = ['Input file: ',filename,'; modified (proxy for created) - ',input_file_date];
options_list{2} = ['Analysed by: ',analysis_file.name,'; ',Rev, '; modified - ',analysis_file.date];
options_list{3} = 'Options: ';
options_names = fieldnames(header_vairables_to_write);
for x = 1 : length(options_names)
    options_list{3} = [options_list{3}, options_names{x}, ' = '];
    if isnumeric(header_vairables_to_write.(options_names{x})) == 1
        options_list{3} = [options_list{3}, num2str(header_vairables_to_write.(options_names{x}))];
    else
        options_list{3} = [options_list{3}, header_vairables_to_write.(options_names{x})];
    end
    if x <= length(options_names) - 1;
        options_list{3} = [options_list{3}, '; '];
    end
end
options_list{4} = 'Sine_fit file version: 2.0'; %this is the output file version number, not the version number of the writing script.


%% Polynomial surface coefficients or background

if nargin > 3

    BGtype = varargin{1};
    Coefs  = varargin{2}; 
%     if nargin > 5
%         Coefs_err  = varargin{3};
%         err_write = 1;
%     else 
%         err_write = 0;
%     end
    
    switch BGtype
        case 'poly surf'
            
            BG_list{1} = 'Polynomial surface coefficients';
            num_str = [7,4];
            
            for x = 1: size(Coefs,1)
                
                Coef_list = Coefs(x,:);
                write_format = [];
                for y = 1:length(Coef_list)
                    %force output string to always be the same length and change
                    %number of decial points to suit
%                     write_form = WriteForm(Coef_list(y), num_str, 2);
                    write_form = ['%', num2str(num_str(1)), '.', num2str(num_str(2)), 'e'];
%                     write_form = ['%10.8e'];
                    write_format = [write_format, write_form];
                    if y ~= length(Coef_list)
                        write_format = [write_format, ', '];
                    end
                        
                    
                end
                BG_list{end+1} = sprintf(['Box%i :  [', write_format,']'], x, Coef_list(:));
            
%                 if err_write == 1
%                     errs_list = Coefs_err(x,:);
%                     write_format = [];
%                     for y = 1:length(errs_list)
%                         %force output string to always be the same length and change
%                         %number of decial points to suit
%                         write_form = WriteForm(errs_list(y), num_str, 2);
%                         write_format = [write_format, ', ', write_form];
%                         
%                     end
%                     BG_list{end+1} = sprintf(['Box%i err:  [', write_format,']'], x, errs_list(:));
%                 end
            end
            
        otherwise
    end
    
    BG_write = 1;
else
    BG_write = 0;
end


%% make output data array
finalised = 0;
max_length = [25, 8];
out = {};
written = zeros(size(variables_to_write,1),1);

x = 1;
fields = fieldnames(variables_to_write);
while finalised == 0
    
    %output data line 1: box/pair numbers:
    if x == 1
        if sum(strcmpi(fields, 'number_boxes')) == 1
            write_field = 'number_boxes';
            write = 1;
        elseif sum(strcmpi(fields, 'number_pairs')) == 1
            write_field = 'number_pairs';
            write = 1;
        else
            error('Write script cannot find the number of pairs/boxes')
        end
        if strcmpi(header_vairables_to_write.length_type, 'pair') == 1 || strcmpi(header_vairables_to_write.length_type, 'anelastic') == 1 
            lead = {['Pair', ','], ['-', ',']};
        elseif strcmpi(header_vairables_to_write.length_type, 'single') == 1
            lead = {['Box number', ','], ['-', ',']};
        else
            error('WriteSineFit cannot determine how the boxes are combined')
        end
        number_columns = variables_to_write.(write_field){:};
        num_write = 1:number_columns;
        write_type = 1;
        
    %output data line 2: horizonal position:    
    elseif x == 2
        loc = strncmpi(fields, 'horizontal', 10);
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Horizontal position,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            error('Write script cannot find the horizontal position in image')            
        end

    %output data line 3: period.
    elseif x == 3
        loc = strcmpi(fields, 'period');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Period,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            if length(num_write)==1
                num_write = repmat(num_write,1,number_columns);
            end
            write_type = 2;
        else
            error('Write script cannot find the period')            
        end
        
    %output data line: period error (if it exists).
    elseif x == 4        
        loc = strcmpi(fields, 'period_err');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Period error,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            if length(num_write)==1
                num_write = repmat(num_write,1,number_columns);
            end
            write_type = 2;
        else 
            write = 0;
        end

    %output data line: phase.
    elseif x == 5
        loc = strcmpi(fields, 'phases');
        if any(loc)==0
            loc = strcmpi(fields, 'phase');
        end
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Phase,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            error('Write script cannot find the Phase')            
        end
        
    %output data line: phase error (if it exists).
    elseif x == 6
        loc = strcmpi(fields, 'phases_err');
        if any(loc)==0
            loc = strcmpi(fields, 'phase_err');
        end
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Phase error,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            write = 0;        
        end
                
    %output data line: amplitude.    
    elseif x == 7
        loc = strcmpi(fields, 'amplitudes');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Amplitude,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            error('Write script cannot find the amplitudes')            
        end
    
    %output data line: amplitude error (if it exists).
    elseif x == 8
        loc = strcmpi(fields, 'amplitudes_err');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Amplitude error,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            write = 0;
        end
            
    %output data line: starting length (if exists).
    elseif x == 9
        loc = strncmpi(fields, 'reference', 8);
        if any(loc)
            write_field = fields{loc};
            write = 1;
            if ~isempty(strfind(write_field, 'length'))
                lead = {'Reference length,', [variables_to_write.(write_field){2},',']};
            elseif ~isempty(strfind(write_field, 'reference'))
                lead = {'Reference position,', [variables_to_write.(write_field){2},',']};
            else
                error('Write script cannot find the reference position or length')
            end
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            write = 0;
        end
        
    %output data line: used in period fitting (if exists).
    elseif x == 10
        loc = strcmpi(fields, 'SSD');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Sum Squared Differences,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            write = 0;
        end
        
    %output data line: number data points in fit (if exists).
    elseif x == 11
        loc = strcmpi(fields, 'nos');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Number data used in fit,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 1;
        else
            write = 0;
        end

    %output data line: used in period fitting (if exists).
    elseif x == 12
        loc = strcmpi(fields, 'useable');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Used in period fit?,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 1;
        else
            write = 0;
        end
        
    %output data line: chisq (if exists).
    elseif x == 13
        loc = strcmpi(fields, 'Chisq');
        if any(loc)
            write_field = fields{loc};
            write = 1;
            lead = {'Chi Squared,', [variables_to_write.(write_field){2},',']};
            num_write = variables_to_write.(write_field){1};
            write_type = 2;
        else
            write = 0;
        end
        
    else
        finalised = 1;
        break;
        
    end 
        
    %make output string
    if write == 1 && finalised == 0
%         outstring = [];
        outstring = sprintf(['%-', num2str(max_length(1)),'s %-', num2str(max_length(2)),'s'], lead{1}, lead{2});
        for y = 1 : length(num_write)
            
            num_str = [9, 4];
            
            %force output string to always be the same length and change
            %number of decial points to suit
            write_form = WriteForm(num_write(y), num_str, write_type);
                        
%             if write_type == 1 && num_write(y)<10^(num_str(1)-num_str(2))
%                 write_form = ['%', num2str(num_str(1)), '.0f'];
%             elseif write_type == 1 && num_write(y)>=10^(num_str(1)-num_str(2))
%                 write_form = ['%',num2str(num_str(1)), '.2g'];
%             elseif write_type == 2 && num_write(y)<10^(num_str(1)-1) && num_write(y)>10^(num_str(1)-num_str(2)-1)
%                 write_form = ['%', num2str(num_str(1)), '.', num2str(num_str(1)-floor(log10(num_write(y)))-2), 'f'];
%             elseif write_type == 2 && num_write(y)>=10^(num_str(1)-num_str(2))
%                 write_form = ['%', num2str(num_str(1)), '.', num2str(num_str(1)-5), 'g'];
%             elseif write_type == 2 && num_write(y)<=10^(-num_str(2))
%                 write_form = ['%', num2str(num_str(1)), '.', num2str(num_str(2)-1), 'e'];
%             else
%                 write_form = ['%', num2str(num_str(1)), '.', num2str(num_str(2)), 'f'];
%             end
            outstring = [outstring, sprintf([write_form, ', '], num_write(y))];
        end
        out{end+1} = outstring;
    end
    
    %check everything is written
    written(x) = 1;
    
    %increment counter
    x = x + 1;
    
end


%% write outputfile
outfile = fopen(OutFileName,'w');
for x = 1:length(options_list)
    fprintf(outfile,'%s\n', options_list{x}); %-- adds header of file names and used options to the '...sine_fits.txt' file
end
if BG_write == 1
    fprintf(outfile,'--\n'); %gap
    for x = 1:length(BG_list)
        fprintf(outfile,'%s\n', BG_list{x}); %-- adds background coefficients to the '...sine_fits.txt' file
    end
end
fprintf(outfile,'--\n'); %gap
for x = 1 : length(out)
    fprintf(outfile,'%s \n', out{x}); %-- writes data lines to the '...sine_fits.txt' file
    fprintf(1,'%s \n', out{x}); %-- writes data lines to the matlab command line.
end
fprintf(1,'\n'); %-- writes data lines to the matlab command line.
fclose('all');


function out = maxnamelength(in)

for x = 1: length(in)
    l(x) = numel(in{x,1});
end
out = max(l);


function out = cellstrfind(cell_array, search_term)
loc = zeros(size(cell_array));
for x = 1 : numel(cell_array)
    y_n = strfind(cell_array{x},search_term);
%     y_n = strncmpi(cell_array{x},search_term, length(search_term));
    if y_n == 1
        loc(x) = 1;
    end
end
if sum(loc(:)) > 1  %finds and dumps the error parameter
    loc2 = find(loc==1);
    for x = 1 : length(loc2)
        res = strfind(cell_array{loc2(x)},'err');
        if ~isempty(res)
            loc(loc2(x)) = 0;
        end
    end
end
out = find(loc~=0);


function write_form = WriteForm(num_write, format_limits, varargin)
%force output string to always be the same length and change
%number of decial points to suit

if nargin >=3
    write_type = varargin{1};
else
    write_type = 2;
end

%format limits
% str_length = format_limits(1);
% decimal_places = format_limits(2);

if write_type == 1 && num_write<10^(format_limits(1)-format_limits(2))
    write_form = ['%', num2str(format_limits(1)), '.0f'];
elseif write_type == 1 && num_write>=10^(format_limits(1)-format_limits(2))
    write_form = ['%',num2str(format_limits(1)), '.2g'];
elseif write_type == 2 && num_write<10^(format_limits(1)-1) && num_write>10^(format_limits(1)-format_limits(2)-1)
    write_form = ['%', num2str(format_limits(1)), '.', num2str(format_limits(1)-floor(log10(num_write))-2), 'f'];
elseif write_type == 2 && num_write>=10^(format_limits(1)-format_limits(2))
    write_form = ['%', num2str(format_limits(1)), '.', num2str(format_limits(1)-5), 'g'];
elseif write_type == 2 && num_write<=10^(-format_limits(2))
    write_form = ['%', num2str(format_limits(1)), '.', num2str(format_limits(2)-1), 'e'];
else
    write_form = ['%', num2str(format_limits(1)), '.', num2str(format_limits(2)), 'f'];
end


%versions
% v3.0 - October 2017
%   - Added code to save the polynomial surfaces from PhasesSSD.
% v2.1 - June 2016
%   -   Added lines so that if the calling function or its version number
%   has dev/DEV in then dev is added to the outfile name.
% v2.0 - November 2014
%   -   The original version of the file. 

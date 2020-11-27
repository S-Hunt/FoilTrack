% ReadSineFit. Reads *sine_fit*.text files from ImageAnalysis suite.
%       Reads the data from the sinefit file and returns it as a series of
%       arrays to the calling function. The function accepts any number of
%       files at the same time as input but they have to be of the same
%       size and format. for which there is only limited checking of format
%       type.
%
%       Syntax: 
%       qc   = ReadSineFit('open', 'file_name_string1', 'file_name_string2'...)
%       qc   = ReadSineFit('close')
%       data = ReadSineFit('var')
%           where variable can be: 
%               'Number'
%               'Position'
%               'Period'
%               'Period error'
%               'Phases'
%               'Phases error'
%               'Amplitudes'
%               'Amplitudes error'
%               'Lengths'
%               'SSD'
%               'Num points'
%               'Used'
%               'ChiSq'
%               'headers'
%               'polynomials'
%               units --- yet to be implemented.
%
%       if qc == 1 then the action was performed correctly but if
%           qc == 0 something went wrong.
%
%       See also: WRITESINEFIT, Phases, PhasesSSD

% Simon Hunt, July, August 2014 & October 2017

function out = ReadSineFit(choice, varargin)

persistent ALLDATA HEADERS number_files POLYNOMIALS; 

%get output variable
choice_line = {
    'Number';...
    'Position';...
    'Period';...
    'Period error';...
    'Phases';...
    'Phases error';...
    'Amplitudes';...
    'Amplitudes error';...
    'Lengths';...
    'SSD';...
    'Num points';...
    'Used';...
    'Chisq'};




if strcmpi(choice, 'open')
   %% read and organise input data. 
    if ~isempty(ALLDATA)
        out = 1;
        disp('Data files have already been opened')
        return
    end
    
    %organise list of file names
    if ischar(varargin)
        number_files = length(varargin);
        file_names = char(varargin);
    else
        number_files = nargin-1;
        file_names = char(varargin);
    end
    
    % reads in all data files listed and does some initial processing of the data.
    for x = 1:number_files
        in = strtrim(file_names(x,:));
        [all_data(x),~,headers(x)] = importdata(in,',');
        all_numbers(x,:,:) = all_data(x).data;
        text_cont(:,:,x) = all_data(x).textdata;
    end
    
    vers = char(text_cont(4,1,:));
    [~,d] = find(vers == ':');
    versions = cellstr(strtrim(vers(:, d+1:end)));
    num_ver = unique(versions);
    
    if length(num_ver) ~= 1 || length(unique(text_cont(4,1,:))) ~= 1
        error('Input files are of different types')
    end
    
    if strcmpi(vers(1,:), '--') == 1
        
        %reads list of headers
            HEADERS = all_data(1).textdata;
        %reads number of data points
            ALLDATA(:,1,:) = all_numbers(:,1,:);
        %reads radii from all_numbers
            ALLDATA(:,2,:) = all_numbers(:,2,:);
        %reads periods from all_numbers
            ALLDATA(:,3,:) = all_numbers(:,3,:);
            ALLDATA(:,4,:) = NaN;
        %reads phases from all_numbers
            ALLDATA(:,5,:) = all_numbers(:,4,:);
        %reads amplidues from all_numbers
            ALLDATA(:,7,:) = all_numbers(:,5,:);
        %reads line to see if used in period fitting
            ALLDATA(:,12,:) = all_numbers(:,7,:);
        %reads SSD from all_numbers
            ALLDATA(:,10,:) = all_numbers(:,8,:);
        %reads lengths from all_numbers
            ALLDATA(:,9,:) = all_numbers(:,9,:);
        %reads errors on phases and amplitudes (if they exist); reads chi
        %squared and number of data points used
        if size(all_numbers,2) >= 10
            ALLDATA(:,6,:) = all_numbers(:,10,:);
            ALLDATA(:,8,:) = all_numbers(:,11,:);
            ALLDATA(:,11,:) = all_numbers(:,12,:);
            ALLDATA(:,13,:) = all_numbers(:,13,:);
        else
            ALLDATA(:,6,:) = NaN;
            ALLDATA(:,8,:) = NaN;
            ALLDATA(:,11,:) = all_numbers(:,10,:);
            ALLDATA(:,13,:) = all_numbers(:,11,:);
            warndlg('There are no errors for the calculated phases or amplitudes.');
        end        
        
    elseif strcmpi(versions, '2.0') | str2double(versions) >= 2
        
        %reads list of headers
            HEADERS = squeeze(text_cont(1:4,1,:));
            
        %read the polynomials
            %version 2.1 of the file and above may have two sets of '--' to
            %differentiate the parts of the data. the polynomials are
            %between them.
            se = find(ismember(text_cont(:,1), '--')~=0);
            if numel(se) == 2     
                POLYNOMIALS = permute(squeeze(text_cont(se(1)+2:se(2)-1,1,:)),[2,1]);
            end
            for x = 1:numel(POLYNOMIALS)
                st = find(POLYNOMIALS{x} == '[');
                POLYNOMIALS{x} = sscanf(POLYNOMIALS{x}(st+1:end-1), '%f,');
            end
        %reads number of data points
            ALLDATA(:,1,:) = all_numbers(:,1,:);
        
        %reads horizontal position from all_numbers
            row = cellstrfind(text_cont(headers+1:end,1), 'Horizontal');
            ALLDATA(:,2,:) = all_numbers(:,row,:);
        
        %reads periods from all_numbers
            row = cellstrfind(text_cont(headers+1:end,1), 'Period');
            ALLDATA(:,3,:) = all_numbers(:,row,:);
        
            if strcmpi(text_cont(headers+row+1,1), 'Period error') == 1
                ALLDATA(:,4,:) = all_numbers(:,row+1,:);
            else ALLDATA(:,4,:) = NaN;
            end
        
        %reads phases from all_numbers
            row = cellstrfind(text_cont(headers+1:end,1), 'Phase');
            ALLDATA(:,5,:) = all_numbers(:,row,:);
        
            if strcmpi(text_cont(headers+row+1,1), 'Phase error') == 1
                ALLDATA(:,6,:) = all_numbers(:,row+1,:);
            else ALLDATA(:,6,:) = NaN;
            end
        
        %reads amplidues from all_numbers
            row = cellstrfind(text_cont(headers+1:end,1), 'Ampl');
            ALLDATA(:,7,:) = all_numbers(:,row,:);
        
            if strcmpi(text_cont(headers+row+1,1), 'Amplitude error') == 1
                ALLDATA(:,8,:) = all_numbers(:,row+1,:);
            else ALLDATA(:,8,:) = NaN;
            end
        
        %reads lengths from all_numbers
            row = cellstrfind(text_cont(headers+1:end,1), 'Starting length')   ;
            if isempty(row)
                row = cellstrfind(text_cont(headers+1:end,1), 'Reference')   ;
            end
            ALLDATA(:,9,:) = all_numbers(:,row,:);
        
        %reads SSD from all_numbers
            row = cellstrfind(text_cont(headers+1:end,1), 'Sum Squared');
            if ~isempty(row)
                ALLDATA(:,10,:) = all_numbers(:,row,:);
            end
                
        %reads number of data points in fit
            row = cellstrfind(text_cont(headers+1:end,1), 'Number data')   ;
            if ~isempty(row)
                ALLDATA(:,11,:) = all_numbers(:,row,:);
            end
            
        %reads values showing if used in phase fitting
            row = cellstrfind(text_cont(headers+1:end,1), 'Used')   ;
            if ~isempty(row)
                ALLDATA(:,12,:) = all_numbers(:,row,:);
            end
        
        %reads Chi Squared values from all_numbers
            row = cellstrfind(text_cont(headers+1:end,1), 'Chi');
            if ~isempty(row)
                ALLDATA(:,13,:) = all_numbers(:,row,:);
            end
        
    else
        error('File version is not recognised')
    end

    out = 1; 


elseif strcmpi(choice, 'close')
    %% close data
    clear ALLDATA HEADERS POLYNOMIALS
    
    out = 1; 


elseif strcmpi(choice, 'headers')
    
    out = HEADERS;
    
elseif strcmpi(choice, 'polynomials')
    
    out = POLYNOMIALS;
    
elseif sum(strcmpi(choice, choice_line)) >= 1
    %% get output data
    if isempty(ALLDATA)
        error('The data files have not been opened')
    end
    
    str_compare = strcmpi(choice, choice_line);
    location = find(str_compare == 1);
    
    if isempty(location)
        error('Unknown data requested -- need to change choice switch')
    else
        out(:,:) = ALLDATA(:,location,:);
        
        %sorts output data
        if isnan(out) %if NaN then makes array a single value
%             out = NaN;
        elseif strcmpi(choice, choice_line{1}) %if number then want single number for boxes/pairs
            out = out(end);
        elseif number_files == 1 %reorientates data if there is only one file - otherwise they are the wrong way round
            out = out';
        end
    end
    
else
    
    out = 0;
    
end


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


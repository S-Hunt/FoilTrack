% ReadPositionChangeFile.
%   Read position change files made by ImageAnalysis
%
%   syntax: [position_change, time_stamps, box_positions, number_boxes, ...
%            number_images] = ReadPositionChangeFile('file_name')
%   
%           headers = ReadPositionChangeFile('file_name', 'header');
%
%   See Also: ImageAnalysis, Phases, PhasesSSD

% Simon Hunt 2014, 2017-2018

function varargout = ReadPositionChangeFile(file_name, varargin)


def = 'positions';
pass = {'warn off'};

% prase varargin...
iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case 'header'
            def = varargin{iarg};
            iarg = iarg + 1;
        case 'positions'
            def = varargin{iarg};
            iarg = iarg + 1;
        case 'warn off'
            pass = [pass, 'warn off'];
            iarg = iarg + 1;
        otherwise
            error(['Unknown option: ' varargin{iarg}]);
    end
end

[data,delimiterOut,headerlinesOut] = importdata(file_name);
if strcmpi(delimiterOut, '	') && headerlinesOut == 1
    %output from multifit.
    
    if strcmpi(def, 'headers')
        error('This will not work! Give up and go home')
    end
    time_stamps = data.data(:,1);
    num_peaks = (size(data.data,2)-3)/2;
    peak_cols = 4:2:(size(data.data,2));
    foil_change_data = data.data(:,peak_cols);
    

    box_positions = zeros(4,num_peaks);
    box_positions(1,:) = mean(foil_change_data,1);
    box_positions(2,:) = box_positions(1,:);
    varargout = {foil_change_data, time_stamps, box_positions, num_peaks, ...
        size(foil_change_data,2)};
    %keyboard
else
    
    % reads in the csv output from the image_analysis.m script
    data = importdata(file_name, ',', 5);
    
    
    
    
    switch def
        case 'positions'
            
            %determine how many header rows there are in data.textdata
            %this is needed because of a change in how importdata works between
            %matlab 2012 and 2017b -- it removes empty lines from text data.
            all_numbers = data.data;
            size_all_numbers = size(all_numbers);
            num_ims = size_all_numbers(1)-4;
            
            % extracts the time stamps of the images from the input file.
            time_stamps = data.textdata(:,3);
            %time_stamps = time_stamps(10:length(time_stamps)); %old line now defunct
            time_stamps = time_stamps(size(time_stamps,1)-num_ims+1:end); %replacement line
            time_stamps = char(time_stamps);
            time_stamps = str2num(time_stamps); %#ok<ST2NM> - this is an array
            
            
            % sorts the numeric data into the box positions in the
            % image and their displacements with time.
            foil_change_data = all_numbers(5:size_all_numbers(1), ...
                1:size_all_numbers(2));
            box_positions = all_numbers(1:4, 1:size_all_numbers(2));
            
            % box_positions(1,N) is the top of the Nth box, (2,n) is the
            % bottom, (3,n) is the left, (4,n) is the right.
            
            % define variables describing amount of data
            number_images = length(time_stamps);
            number_boxes = size(foil_change_data, 2);
            
            %get image numbers.
            for x = 9:size(data.textdata,1)
                first_im{x-8,:} = FileTitleInformation(data.textdata{x,1}, pass{:});
                second_im{x-8,:} = FileTitleInformation(data.textdata{x,2}, pass{:});
                im1_num(x-8) = first_im{x-8}{5};
                im2_num(x-8) = second_im{x-8}{5};
            end
            im_nums = [im1_num', im2_num'];
            
            %output
            varargout = {foil_change_data, time_stamps, box_positions, number_boxes, ...
                im_nums};
                %number_images};
        case 'header'
            
            %parse the headers
            
            %inputfile
            sep = find(data.textdata{1}==':',1,'first');
            str = data.textdata{1}(sep+1:end);
            head.file = strtrim(str);
            
            %Analysis code.
            sep = find(data.textdata{2}==':');
            sep2 = find(data.textdata{2}=='.');
            
            head.analysis.script{1} =  strtrim(data.textdata{2}(sep(1)+1:sep2(1)+2)); %after first colon to after first fullstop + 2
            head.analysis.script{2} =   strtrim(data.textdata{2}(strfind(data.textdata{2},'and')+4:sep(end)-8));% after 'and' to 8 before last colon
            
            head.analysis.modified{1} =  strtrim(data.textdata{2}(find(data.textdata{2}=='(')+9:find(data.textdata{2}==')')-1)); % after '('+9 to before ')'
            
            head.analysis.version{1} =  strtrim(data.textdata{2}(sep(2)+1:find(data.textdata{2}=='(')-1)); %after second colon to before first '('
            head.analysis.version{2} =  strtrim(data.textdata{2}(sep(end)+1:end));% last colon to the end.
            
            
            %analysis options in file.
            sep = find(data.textdata{3}==':');
            sep2 = find(data.textdata{3}==';');
            sep3 = find(data.textdata{3}=='-');
            
            for x = 1:length(sep3)
                
                if x == length(sep3)
                    bigger = length(data.textdata{3}) + 1;
                else
                    bigger = sep2(find(sep2>sep3(x), 1, 'first'))
                end
                arr = [sep, sep2]
                smaller = arr(find(arr < sep3(x), 1, 'last'))
                %             smaller =
                
                
                
                
                f = (data.textdata{3}(smaller+1: sep3(x)-1));
                f = f(~isspace(f));
                f = regexprep(f, '/', '');
                f = regexprep(f, '+', '');
                
                v = strtrim(data.textdata{3}(sep3(x)+1:bigger-1));
                
                head.options.(f) = v;%version{1}
                
            end
            
            
            
            %output
            varargout = {head};
            
            
    end
    
end
    
    
end


%Changes
% v2.1 - June 2018
%   - bug fix the account for differences between Matlab 2012 and 2017b.
% v2 - Dec 2017
%  -parse the header in the position change file.
%v1 - 2014
% - read position change text file but not the header.
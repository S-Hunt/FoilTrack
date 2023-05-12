% WritePositionChange
%   Writes the position changes calcualted by ImageAnalysis. 
%
%   Syntax: 
%       status = WritePositionChange(OutFileName, 'open')
%       status = WritePositionChange(OutFileName, 'header', header_structure)
%       status = WritePositionChange(OutFileName, 'values', ...
%                                              offsets, out_ref, out_nos, out_times )
%       status = WritePositionChange(OutFileName, 'close')
%
%   The file has the following format:
%       first line:   String listing data file name
%       second line:  String listing files performing the analysis
%       third line:   String listing the used Options, which are separated by ';'
%       fourth line:  -- blank
%       lines 5-9:    Position of the boxes (comma separated and lined up with displacements below them)
%       line 10:      -- blank
%       lines 11+:    data formatted as:
%           RefImageName, ComparitorImageName, TimeStamp1, DisplacementOfBox1atTimeStamp1, DisplacementOfBox2atTimeStamp1...
%           RefImageName, ComparitorImageName, TimeStamp2, DisplacementOfBox1atTimeStamp2, DisplacementOfBox2atTimeStamp2...
%           ...
%       end of data:  -- blank line
%
%    See also: ImageAnalysis, READpositionChange

% Simon Hunt May 2016
% $ version 0.1.1 $ May 2016 $
% Code stripped from ImageAnalysis 3.9.1 and edited. 

function out = WritePositionChange(OutFileName, process, varargin)

persistent outfile_handle

switch process
    case {'Open'} %open the file for writing.
        outfile_handle = fopen(OutFileName,'a');
        out = 1; 
        
    case {'OpenNew'} %open the file for writing.
        outfile_handle = fopen(OutFileName,'w');
        out = 1; 
        
    case 'header' %write the header to the displacements file. 
        
        disp('Creating header of the position_change file.')
        
        to_write = varargin{1};
        
        %Header:
        % first line:   Input file: Data_file_name
        % second line:  Analysed by: ImageAnalysis.m version: 1.9 (modified 02-May-2016 15:23:50) and displacement_analysis version: 3.9 
        % third line:   Options: 
        % fourth line:  -- blank
        % lines 5-9:    Position of the  

        % String for fist line of header file -- 
        options_list{1} = ['Input file: ',to_write.file_name];
        
        % Information for second line of header. 
        file_date = dir(to_write.caller);
        
        %reads version of the caller (which is usually ImageAnalysis.m)
        [Ver,content] = FileVersion(to_write.caller);

        %finds which file is used for the displacement and reads is version number.
        subcall = CalledFile(content, '%<-disp analysis row->%');
        Ver2 = FileVersion(which(subcall));
        
        % String for second line of header file -- 
        options_list{2} = ['Analysed by: ',file_date.name,' ',Ver, ' (modified ',file_date.date,') and ', subcall, '.m ', Ver2];
        
        % String for third line of header file --
        %check rotation is present if not make = 0
        if isfield(to_write.boxes, 'Image_rotation')==0
            to_write.boxes.Image_rotation =0;
        end
        options_list{3} = ['Options: Displacement type - ',to_write.variables.analysis_type,...
            '; reference image - #',num2str(to_write.refid),...
            '; image rotation - ',num2str(to_write.boxes.Image_rotation),...
            '; minimum of SSD mode - ',to_write.variables.min_type,...
            '; search +/-',num2str(to_write.variables.search),...
            '; bright spot removal - ', to_write.variables.spot_removal,...
            '; NaN area outside foil - ', to_write.variables.NaN_bg];
        if isfield(to_write.variables ,'begin')
            options_list{3} = [options_list{3},...
                '; begin - ',num2str(to_write.variables.begin),...
                '; cease - ',num2str(to_write.variables.cease),...
                '; window - ',num2str(to_write.variables.window)];
        end
        if isfield(to_write.variables ,'surface_coef')
            options_list{3} = [options_list{3},...
                '; surface coefficients - [',num2str(to_write.variables.surface_coef(1)), ',',num2str(to_write.variables.surface_coef(2)),']'];
        end
        
        
        %write first three lines of the header file.
        for x = 1:length(options_list)
            fprintf(outfile_handle,'%s \n', options_list{x}); %-- adds header of file names and used options
        end
        
        %header line 4. Empty line
        fprintf(outfile_handle,'\n'); %gap
        
        % Write the positions of the boxes to the header
        name1 = 'First_image';
        name2 = 'Second_Image';
        name3 = 'box';
        name5 = 'time-of-second-image';
        
        %introduction line for box positions.
        fprintf(outfile_handle, '%s, %s, %s', name1, name2, name5);
        for x = 1 : length(to_write.boxes.boxX)
            fprintf(outfile_handle, ', %s%1.0f', name3, x);
        end
        fprintf(outfile_handle, ' \n');
        
        %write box positions themselves
        name1 = strcat(to_write.run_name, ', -');
        for times = 1 : 4
            if times == 1
                name2 = 'top   ';
                position = to_write.boxes.boxY(:,1);
            elseif times == 2
                name2 = 'bottom';
                position = to_write.boxes.boxY(:,2);
            elseif times == 3
                name2 = 'left  ';
                position = to_write.boxes.boxX(:,1);
            else
                name2 = 'right ';
                position = to_write.boxes.boxX(:,2);
            end
            
            fprintf(outfile_handle, '%s, %s', name1, name2);
            for x = 1 : size(to_write.boxes.boxX,1)
                fprintf(outfile_handle, ', %8.5f', position(x) );
            end
            fprintf(outfile_handle, ' \n');
        end
        fprintf(outfile_handle, '\n');
        
        
        out = 1;
        

    case 'profile header' %write the header to the displacements file. 
        
        disp('Creating header of the ''profile'' file.')
        
        to_write = varargin{1};
        
        %Header:
        % first line:   Input file: Data_file_name
        % second line:  Analysed by: ImageAnalysis.m version: 1.9 (modified 02-May-2016 15:23:50) and displacement_analysis version: 3.9 
        % third line:   Options: 
        % fourth line:  -- blank
        % lines 5-9:    Position of the  

        % String for fist line of header file -- 
        options_list{1} = ['Input file: ',to_write.file_name];
        
        % Information for second line of header. 
        file_date = dir(to_write.caller);
        
        %reads version of the caller (which is usually ImageAnalysis.m)
        [Ver,content] = FileVersion(to_write.caller);

        %finds which file is used for the displacement and reads is version number.
        subcall = CalledFile(content, '%<-disp analysis row->%');
        Ver2 = FileVersion(which(subcall));
        
        % String for second line of header file -- 
        options_list{2} = ['Analysed by: ',file_date.name,' ',Ver, ' (modified ',file_date.date,')'];
        
        
        %write first two lines of the header file.
        for x = 1:length(options_list)
            fprintf(outfile_handle,'%s \n', options_list{x}); %-- adds header of file names and used options
        end
        
        %header line 3. Empty line
        fprintf(outfile_handle,'\n'); %gap
        
        %introduction line for box positions.
        fprintf(outfile_handle, 'Band averaged over:, X_min,        %i' , to_write.Xmin);
        fprintf(outfile_handle, ' \n');
        fprintf(outfile_handle, '                     X_max,        %i' , to_write.Xmax);
        fprintf(outfile_handle, ' \n');
        
        %header line 6. Empty line
        fprintf(outfile_handle,'\n'); %gap
        
        
        fprintf(outfile_handle, 'File,                     , Time Stamp');
        fprintf(outfile_handle, ' \n');
        
        
        out = 1;
        
        
    case 'values'
        
        if nargin == 6 % position change file
            offsets    = varargin{1};
            out_ref    = varargin{2};
            out_nos    = varargin{3};
            out_times  = varargin{4};
        elseif nargin == 5 %% profile file
            offsets    = varargin{1};
            out_ref    = varargin{2};
            out_nos    = varargin{2}; %make an array so that it has something to pass around. This is a very poor way of doing it but hey-ho.
            out_times  = varargin{3};
        end
            
        num_boxes = size(offsets,2);
        num_lines = numel(out_times);
        
        %make the format string the right size for the number of boxes
        format = ['%s, %s, %5.3f', repmat(', %4.7f',1,num_boxes),' \n'];
        
        for x = 1 : num_lines
            %N.B. out_ref and out_nos need to be cell arrays inorder to cope with any 
            % change in length of the strings.
            
            %if the images are files (i.e. not an netcdf file) cut the directory name from the strings.
            if ~isempty(strfind(out_ref{x},'\')) || ~isempty(strfind(out_ref{x},'/'))
                try
                    [~,a1,b1] = fileparts(cell2mat(out_ref{x}));
                catch
                    [~,a1,b1] = fileparts(out_ref{x});
                end
                ref_name = [a1,b1];
                try
                    [~,a2,b2] = fileparts(cell2mat(out_nos{x}));
                catch
                    [~,a2,b2] = fileparts((out_nos{x}));
                end
                comp_name = [a2,b2];
            elseif isnumeric(out_ref{x}) == 1
%                 ref_name = char(out_ref{x,:});
%                 comp_name = char(out_nos{x,:});
                ref_name = sprintf('%d', out_ref{x});
                comp_name = sprintf('%d', out_nos{x});
            else
                ref_name = char(out_ref{x});
                comp_name = char(out_nos{x});                
            end
            if nargin == 6 
                fprintf(outfile_handle, format, ref_name, comp_name, out_times(x), offsets(x,:));
            elseif nargin == 5
                fprintf(outfile_handle, format, ref_name, '-', out_times(x), offsets(x,:));
            end
        end
        
        out = 1;
        
        
    case 'close'
        
        %write empty lines so if start data 
        fprintf(outfile_handle, ' \n');
        
        out = fclose(outfile_handle);
        if out == 0
            out = 1;
        end
        
end %switch process

%error out if the process does not work
if out ~= 1
    error(['The WriteFunction has not worked for ''',process,''''])
end


end %WritePositionChange


function [version, content] = FileVersion(file)

%reads version of this file.
fid = fopen(file);
content = fscanf(fid,'%c');
dollars = find(content == '$');
version =  strtrim(content(dollars(1)+1:dollars(2)-1));
fclose(fid);

end %FileVersion


function disp_call_name = CalledFile(content, search_string)

where = strfind(content, search_string);
% where = strfind(content, '%<-disp analysis row->%');
line_content = content(where(2)-200:where(2));
line_begin = find(line_content == '=');
line_end = find(line_content == '(');
disp_call_name = line_content(line_begin(end)+2:line_end(end)-1); %file name of displacement analysis file
        
end % CalledFile


%%versions
% v0.1.1 May 2017
%   -   Added catch for image rotation to make 0 if it is not present in
%   the variables.
%v0.1 May 2016
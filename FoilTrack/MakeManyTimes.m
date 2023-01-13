% MAKEMANYTIMES.
% Generates the batch files for bulk processing of image data from synchrotron experiments.
% The batch files call variously: ImageAnalysis, Phases, PhasesSSD, KappaSolve 
% Syntax:
%       MakeManyTimes('target directory')
%       MakeManyTimes('target directory', 'experiment type',...)
%       MakeManyTimes('target directory', 'order',...)
%       MakeManyTimes('target directory', 'image format',...)
%       MakeManyTimes('target directory', 'proc',...)
%
%   - Possible experiemnt types are 'rheology', 'anelasticity' and 'td'. If there is no experiment type 
%        assumes rheology experiment and only makes manytimes.m
%   - Possible orders are 'acsend' and 'descend'.
%   - Sort: followed by options 'load', 'temperature', 'period' in order of presidence
%   - Possible image formats are 'tiff', 'tif', 'bmp', 'netcdf', 'nc' and 'lst'. If none are 
%        present 'netcdf' is assumed.
%   - Set procedure to be undertaken in manytimes.m. ...'proc', 'setting'.
%        For the list of options see ImageAnalysis.m
%
%   See Also MakeManyLstFiles, MakeSingleLstFile, ImageAnalysis

%   Simon Hunt 2007 - 2019
%   Version 2.2.1  -- May 2019

function MakeManyTimes(target_dir, varargin)

%process varargin and defaults
direction = 'ascend';
image_type = '';
expt_type = '';
sort_by = {};
open_to_edit = 0;
SSD = 0;
proc = 'boxes';

iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case {'descend', 'ascend'}
            direction = varargin{iarg};
            iarg = iarg + 1;
        case {'rheology', 'anelasticity', 'td'}
            expt_type = varargin{iarg};
            iarg = iarg + 1;
        case {'tiff', 'tif', 'bmp', 'netcdf', 'nc', 'lst'}
            image_type = varargin{iarg};
            iarg = iarg + 1;
        case {'edit'}
            open_to_edit = 1;
            iarg = iarg + 1;
        case {'ssd'}
            SSD = 1;
            iarg = iarg + 1;
        case 'proc'
            proc = varargin{iarg+1};
            iarg = iarg + 2;
        case {'sort'}
            va = {'load', 'temperature', 'period', 'timestamp'};
            %this should be a while loop but the while loop would not work
            %when I wrote it. Hence it is a for loop.
            for x = 1:3
                if iarg+x <= length(varargin) && ismember(varargin{iarg+x}, va)
                    sort_by = [sort_by; varargin{iarg+x}];
                end
            end
            iarg = iarg + 1 + length(sort_by);
        otherwise
            error(['Unknown option: ' varargin{iarg}]);
    end
end

%check trget directory exists and is correctly formatted.
if exist(target_dir, 'dir') ~= 7
    error('The target folder does not exist');
end

if target_dir(end) ~= filesep
    target_dir = strcat(target_dir,filesep);
end

%% list image files

%determine what to list in the batch files.
%determine if files are netcdfs or single image files.
netcdfs = sum(strcmpi(image_type, {'netcdf', 'nc'}));
% lsts = strcmpi(image_type, 'lst');

if netcdfs == 1 %the file type has been set to be netcdf
    file_search = '*.nc';
    dir_search = target_dir;
    
else %the file type has not been set or is set not to be a netcdf.
    
    %first look for netcdf files in the target directory.
    allimages_nc = dir(strcat(target_dir,'*.nc')); %list netcdf files
    allimages_lstlocal = dir([cd,filesep,'*.lst']); %check current directory for lst files (there should not be any)
    allimages_lst = dir(strcat(target_dir,'*.lst')); %check target directory for lst files (there should not be any)
    if ~isempty(allimages_nc)
        file_search = '*.nc';
        dir_search = target_dir;
    elseif ~isempty(allimages_lstlocal)          
        file_search = '*.lst';
        dir_search = [cd,filesep];
    elseif ~isempty(allimages_lst)      
        file_search = '*.lst';
        dir_search = target_dir; 
    else
        
        fprintf('No netcdf or list files have been found.\n Make the .lst files for the target images.\n')
        
        %arguemnts to pass to MakeManyLstFiles
        args_to_pass = {direction};
        if ~isempty(image_type) && strcmpi(image_type, 'lst') == 0
            args_to_pass = {args_to_pass{:}, image_type};
        end
        
        %make lst file in current directory
        MakeManyLstFiles(target_dir, args_to_pass{:});
%         
%         %search for the listings again
%         allimages = dir([dir_search,'.lst']);
%         
        file_search = '*.lst';
        dir_search = [cd,filesep]; 
    end

end
    
%list all required files.
allfiles = dir(strcat(dir_search,file_search));

%if allimages is empty 
if isempty(allfiles)
    error('No recognised image or .lst files have been found.') 
end
% if isempty(allimages) && netcdfs == 0
%     fprintf('List files have not been found. \n Make the .lst files for the target images.\n')
%     
%     %arguemnts to pass to MakeManyLstFiles
%     args_to_pass = {direction};
%     if ~isempty(image_type) && strcmpi(image_type, 'lst') == 0
%         args_to_pass = {args_to_pass{:}, image_type};
%     end
%     
%     %make lst file in current directory
%     MakeManyLstFiles(target_dir); 
% 
%     %search for the listings again
%     allimages = dir([dir_search,'.lst']);
%     
% elseif isempty(allimages) && netcdfs == 1    
%     error(['No netcdf image files were found in ', dir_search]) 
% end

num_files = length(allfiles);

fprintf('Found %i %s files to add to batch files.\n', num_files, file_search(3:end))

%order the netcdf files according to the incremental counter in the names 
for x = 1:num_files
    name = allfiles(x).name;
    titles(x,:)  = FileTitleInformation(name, 'warn off');
    nos(x) = titles{x,5};
%     root_length = find(name == '.');
%     order_num = name(1,root_length(end)-3:root_length(end)-1);
%     try 
%         nos(x) = str2num(order_num);
%     catch
%         nos(x) = NaN;
%     end
end

% sorts the data by collection or forced order.
if strcmpi(sort_by, 'timestamp') == 1 
    %sort by the timestamps of the first image in the series
    for x = 1:num_files
        if netcdfs == 1 %if netcdf
            error('The option to sort netcdf files by timestamp has not been created')
            % FIXME
        else %if not netcdf (assume lst file)
            fid = fopen(allfiles(x).name);
            temp = fgetl(fid);
            coms = find(temp == ',');
            C(x) = str2num(strtrim(temp(coms(1)+1:coms(2)-1)));
% C(x) = textscan(fid, '%*s %n,', 1);
            fclose(fid);
        end
    end
%     C = cell2mat(C);
    [~,order] = sort(C);
elseif ~isempty(sort_by)
    %sort the data by the order in sort_by
    
    %replace strings by colum numbers.
    sort_by = strrep(sort_by', 'load', '2');
    sort_by = strrep(sort_by, 'temperature', '3');
    sort_by = strrep(sort_by, 'period', '4');
    sort_by = str2num(cell2mat(sort_by'));
    
    %make array of order and values. 
    ltp = [(1:num_files)', cell2mat(titles(:,2:4))];
    for x = length(sort_by):-1:1
        ltp = sortrows(ltp, sort_by(x));
    end
    order = ltp(:,1);
else
    %sort by the counting number of the data sets (added to NetCDF files by default). 
    order_list(:,1) = 1:num_files;
    order_list(:,2) = nos;
    order_sort = sortrows(order_list,2);
    order = order_sort(:,1);
end

% order_list(:,2) = 1:num_files;
% order_sort = sortrows(order_list);
if strcmpi(direction, 'descend') == 1
%     order = order_sort(:,2);
    order = order(end:-1:1); %reverses the order if descending.
% else % strcmpi(direction, 'ascend') == 1
%     order = order_sort(:,2);    
end
allfiles(1:num_files,:) = allfiles(order,:);


%% make manytimes output file
processname = 'ImageAnalysis';
opening = '( strcat(data_dir,''';
mid = '''), ';
ending = ');';
    
name2 = 'process';
name3 = 'type';
    
x_from = 1;

outfile = fopen('manytimes.m', 'w');

fprintf(outfile, 'process = ''%s''; %% See ImageAnalysis for details of options can use\n', proc);
fprintf(outfile, 'data_dir = ''%s''; \n \n', dir_search);

x_to = num_files;
for x = x_from : x_to
    name1 = allfiles(x).name;

    paste = horzcat(processname, opening, name1, mid, name2, ending);
    fprintf(outfile, '%s \n', paste);
end
    
fprintf(outfile, '\nclose all;');
fclose(outfile);

fprintf('Written manytimes.m to the current directory \n');

%% make manytimes2 output file

if sum(strcmpi(expt_type, {'anelasticity', 'td'})) ~=0
    
    processname = 'Phases';
    opening = '(''';
    ending = ''');';
    file_ending = '_position_change.txt';
    
    outfile = fopen('manytimes2.m', 'w');
    for x = 1 : length(allfiles)
        name = allfiles(x).name;
        root_length = find(name == '.',1,'last');
        outname = strcat(name(1:root_length-1),file_ending);
        
        paste = horzcat(processname, opening, outname, ending);
        fprintf(outfile, '%s \n', paste);
    end
    fprintf(outfile, '\nclose all;');
    fclose(outfile);
    fprintf('Written manytimes2.m to the current directory \n');
    
    if SSD == 1
        processname = 'PhasesSSD';
        opening = '(''';
        ending = ''');';
        file_ending = '_SSD.mat';
        
        outfile  = fopen('manytimes2_SSD.m', 'w');
        for x = 1 : length(allfiles)
            name = allfiles(x).name;
            root_length = find(name == '.',1,'last');
            outname = strcat(name(1:root_length-1),file_ending);
            
            paste = horzcat(processname, opening, outname, ending);
            fprintf(outfile, '%s \n', paste);
        end
        fprintf(outfile, '\nclose all;');
        fclose(outfile);
        
        fprintf('Written manytimes2_SSD.m to the current directory \n');
    end
end

%% make manytimes2parallel output file
if sum(strcmpi(expt_type, {'anelasticity', 'td'})) ~=0
    done = 1;
    opening = '    ''';
    ending = ''';';
    while done ~= 0
        if done == 1
            processname = 'Phases';
            file_ending = '_position_change.txt';
            outfiel_name = 'manytimes2parallel.m';
        elseif done == 2
            processname = 'PhasesSSD';
            file_ending = '_SSD.mat';
            outfiel_name = 'manytimes2parallel_SSD.m';
        end
        outfile = fopen(outfiel_name, 'w');
        fprintf(outfile, '%% Parallel processing of data to find phases.\n');
        fprintf(outfile, '%%All the functionality comes after the list of file names.\n\n');
        fprintf(outfile, 'file_names = {\n');
        for x = 1 : length(allfiles)
            name = allfiles(x).name;
            root_length = find(name == '.',1,'last');
            outname = strcat(name(1:root_length-1),file_ending);
            
            paste = horzcat(opening, outname, ending);
            fprintf(outfile, '%s \n', paste);
        end
        fprintf(outfile, '          };\n\n');
        %fprintf(outfile, 'if matlabpool(''size'') == 0\n');
        %fprintf(outfile, '    matlabpool open;\n');
        %fprintf(outfile, 'end\n');
        fprintf(outfile, 'parpool\n');
        fprintf(outfile, 'loops = length(file_names);\n');
        fprintf(outfile, 'parfor x = 1 : loops\n');
        fprintf(outfile, ['     ',processname,'(file_names{x});\n']);
        fprintf(outfile, 'end\n');
        %fprintf(outfile, 'matlabpool close;\n');
        fprintf(outfile, 'close all;');
        fclose(outfile);
    
    
        fprintf(['Written ',outfiel_name, ' to the current directory \n']);
    
        if done == 1 && SSD == 1
            done = 2;
        elseif done == 2
            done = 0;
        else
            done = 0;
        end            
    end
end

%% make manytimes3 output file
if strcmpi(expt_type, 'td') == 1
    
    processname = 'Kappa_solve';
    opening = '(''';
    midstart = '''';
    midend = ''', ...';
    ending = ''');';
    file_ending = '_sine_fits.txt';
    
    for x = 1:num_files
        name = allfiles(x).name;
        conditions = FileTitleInformation(name);
        force(x) = conditions{2};
        temp(x) = conditions{3};
    end
    
    uforce = unique(force);
    utemp = unique(temp);
    if strcmpi(direction, 'descend') == 1
        uforce = fliplr(uforce);
        utemp = fliplr(utemp);
    end
    f = 1;
    t = 1;
    done = 0;
    
    outfile = fopen('manytimes3.m', 'w');
    while done == 0
        this_time = find(force == (uforce(f)) & temp == utemp(t)); %finds each unique paring of load and temp
        if isempty(this_time) ~= 1
            fprintf(outfile, '%s', [processname, opening]);
            for x = 1 : length(this_time)
                name = allfiles(this_time(x)).name;
                root_length = find(name == '.',1,'last');
                outname = strcat(name(1:root_length-1),file_ending);
                fprintf(outfile, '%s', outname);
                if x == length(this_time)
                    fprintf(outfile, '%s\n\n', ending);
                else
                    fprintf(outfile, '%s\n    %s', midend, midstart);
                end
            end
        end
        % switches to work through all combinations of load and temperature
        if t == length(utemp) && f == length(uforce)
            done = 1;
        elseif t ~= length(utemp)
            t = t + 1;
        elseif t == length(utemp)
            f = f + 1;
            t = 1;
        end
    end
    fprintf(outfile, '\nclose all;');
    fclose(outfile);

    
    fprintf('Written manytimes3.m to the current directory \n');
end


%% open files
if open_to_edit == 1
    edit manytimes 
    if sum(strcmpi(expt_type, {'anelasticity', 'td'})) ~=0
        edit manytimes2 manytimes2parallel
    end
    if strcmpi(expt_type, 'td') == 1
        edit manytimes3
    end
end

clear;

%version chages
% v2.2.1 -- April 2019
%   - chenged matlabpool calls to parpool
% v2.2 -- May 2017
%   - Added option to change the function in manytimes.m
%   - Added option to make batch files to feed SSD outputs into PhaseSSD.
% v2.1 -- Jan 2017
%   - Added option to sort by condition.
%   - Added another/better option to sort by time stamp of first image in the files. 
% V 2.0 - 2016
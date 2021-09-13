% MakeSingleLstFile(directory, opts...).
%
% Makes a *.lst file for the images in the target directory.
% In the list file is the list of images and their associated time stamps.
% Without the switch for 'image format' the script uses the most common
% file type in the directory. 
%
% Syntax:
%       MakeSingleLstFile('target directory')
%       MakeSingleLstFile('target directory', 'image format')
%       MakeSingleLstFile('target directory', 'order')
%       MakeSingleLstFile('target directory', 'temperature')
%       MakeSingleLstFile('target directory', 'logfile', fname)
%       MakeSingleLstFile('target directory', 'replace', fname, [Nkeep])
%       MakeSingleLstFile('target directory', 'plot')
%
%    'image format': recognsed image formats are 'tif', 'tiff', 'edf' and
%                    'bmp'. If none is specified 'tif' is the default.
%                       'nxs' as an image format also requires a second
%                       option naming the camera (for Diamond *.nxs
%                       files). e.g. 'nxs', 'camFloat2'. This is to filter
%                       the XRD or radiography images from each other. 
%    'order':        'ascend' or 'descend' to determine the order the images are
%                    listed in. 'ascend' is the default.
%    'temperature':  Option to read tempeartures from image files and list in the lst file
%    'logfile':      Option to read time stamps from a logfile. Values 0, 1 or 2.
%                           0: - void, does nothing,
%                           1: - logfile [match file names with log file]
%                           2: - replace [match save time with closest time in logfile]
%    'replace':      Option to read time stamps from another file, this just replaces the first 
%                        N time stamps by postion in the list and then minimises the differences 
%                        for the rest.
%                   fname can take the format asdf->fgjk/*.log -- in which
%                   case the script replaces asdf with fghk in the target
%                   directory string and appends everything after the file
%                   separator.
%    'plot':          Plot the time stamps vs. frame number for the data
%                     set.
%
%    For both logfile and replace options a relative location is assumed for the log file 
%       if the string starts with './', '.\', '../' or '..\'.
%
%   See Also MakeManyTimes, MakeManyLstFiles

%	Version 1.4.2 -- 18th Feb 2019
%   Simon Hunt 2016-2019
% This script used to be called makemanytimes_tiff.

function  MakeSingleLstFile(target_dir, varargin)

%process varargin and defaults
plot_on = 0;
direction = 'ascend';
type = 'tif';
add_temp = 0;
add_time = 1;
power_file = [];
logfile = [];
ts = 0; %time stamp replacement settings: 0 = void; 1 = logfile [match file names with log file]; 2 = replace [match save time with closest time in logfile]
rep = [];
camera = 'camFloat2';

iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case {'descend', 'ascend'}
            direction = varargin{iarg};
            iarg = iarg + 1;
        case {'tiff', 'tif', 'bmp', 'edf', 'mat'}
            type = varargin{iarg};
            iarg = iarg + 1;
        case {'nxs'}
            type = varargin{iarg};
            camera = varargin{iarg+1};
            iarg = iarg + 2;
        case {'logfile'}
            logfile = varargin{iarg+1};
            ts = 1;
            iarg = iarg + 2;
        case 'plot'
            plot_on = 1;
            iarg = iarg+1;
        case {'replace'}
            logfile = varargin{iarg+1};
            ts = 2;
            iarg = iarg + 2;
            if numel(varargin)>=iarg && isnumeric(varargin{iarg})
                rep = varargin{iarg};
                iarg = iarg + 1;                
            end
%             error('Logfile input is not implemented.')
        case {'temperature'}
            add_temp = 1;
            iarg = iarg + 1;
        case {'power'}
            power_file = varargin{iarg+1};
            iarg = iarg + 2;
        case {'no time'}
            add_time = 0;
            iarg = iarg + 1;
        otherwise
            error(['Unknown option: ' varargin{iarg}]);
    end
end

%check trget directory exists and is correctly formatted.
if exist(target_dir, 'dir') ~= 7
    error('The target folder does not exist');
end

% make sure target_dir has only the correct file separator in it.
target_dir = strrep(target_dir, '/', filesep);
target_dir = strrep(target_dir, '\', filesep);

% make sure target_dir ends with a file separator.
if target_dir(end) ~= filesep
    target_dir = strcat(target_dir,filesep);
end


% get all images
allimages = dir(strcat(target_dir, '*.', type, '*'));
if isempty(allimages) 
    fprintf(' There are no images with extenstion *.%s. Therefore try using most common file format\n', type)
    
    allimages = dir(target_dir);
    
    %find the most common file type and list only those
    for x = 1:length(allimages)
        if allimages(x).isdir == 1
            ext{x} = 'dir';
        else
            [~,~,ext{x}] = fileparts(allimages(x).name);
        end
    end
    endings = unique(ext);
    
    %remove list files and directorys as possibilities.
    any_dir = ~cellfun(@isempty,strfind(endings,'dir')); %find if there are any dir entries in the list
    any_lst = ~cellfun(@isempty,strfind(endings,'lst')); %find if there are any lst entries in the list
    any_chi = ~cellfun(@isempty,strfind(endings,'chi')); %find if there are any chi entries in the list
    any_db  = ~cellfun(@isempty,strfind(endings,'db'));  %find if there are any db entries in the list
    endings(any_dir) = []; %remove 'dir' as an option
    endings(any_lst) = []; %remove '.lst' as an option
    endings(any_chi) = []; %remove '.lst' as an option
    endings(any_db)  = []; %remove '.db' as an option
    if isempty(endings)
        fprintf(' No files were found in the directory. No *.lst file can be made.\n')
        return
    end
    
    %find most common extension.
    for x = 1:length(endings)
       count(x) =  sum(ismember(ext, endings(x)));
    end
    [~,use] = max(count);
    type = endings{use};
    allimages = allimages(ismember(ext, endings(use)));
end

%strip '.' from type string if it is there
type(type=='.') = [];


%determine filename ending
if strcmpi(type, 'edf') == 1
    ending = 'edflst';
elseif strcmpi(type, 'nxs') == 1
    ending = 'nxslst';
else
    ending = 'lst';
end



%NXS images
%search though the images and determine which are to be listed.
discard = [];
if strcmpi(type,'nxs') == 1
    fprintf(' Searching through %i *.nxs files to determine which are X-radiographs.\n', length(allimages))
    for x = 1:length(allimages)
        cam{x} = Image_Functions_NXSlist('camera', [target_dir, allimages(x).name]);
        
        ind(x) = strcmp(cam{x}, camera);
    end
    allimages(~ind) = [];
end




num_images = length(allimages);

fprintf(' Found %i *.%s files to add to .%s file.\n', num_images, type, ending)

%get file title information and sort the input images by number.
for x = 1:num_images
    file_parts(x,:) = FileTitleInformation(allimages(x).name, 'warn off');
end
[~,ord] = sort(cell2mat(file_parts(:,5)), direction);
file_parts = file_parts(ord,:);
allimages = allimages(ord);

%read logfile (if required)
if ts ~= 0
    
    if ~isempty(strfind(logfile, '->')) %replace a sting in the target dir and append another part of string.
        
        %find replacement indicator '->'
        repl = strfind(logfile, '->');
        
        %find seperator indicator
        repl2 = [strfind(logfile, '/'), strfind(logfile, '\')];
        if isempty(repl2)
            repl2 = length(logfile)
        else
            repl2 = repl2(1);
        end
           
        %get string to add to end
        if repl2+1 < length(logfile)
            add_str = logfile(repl2:end);
        else
            add_str = [];
        end
        
        logfile = [strrep(target_dir, logfile(1:repl-1), logfile(repl+2:repl2-1)), add_str];
        
    elseif strcmpi(logfile(1), '.') == 1 %just append the name to the search
        logfile = [target_dir, logfile];
    end
    [ImTimes ImNames] = ReadLogFile(logfile);
end


%get time stamps of images
for i = 1:num_images
    im_name = [target_dir, allimages(i).name];
    try
        im_info = imfinfo(im_name);
    catch
        im_info = [];
    end
%     if ts == 1 %get the times of the images from the logfile.        
%         loc = find(ismember(ImNames,  allimages(i).name));
%         if isempty(loc)
%             error_str = ['No time stamp has been found for ',allimages(i).name,' in the logfile.'];
%             error(error_str);
%         end
%         im_times(i) = ImTimes(loc);       

%else
    if isfield(im_info, 'UnknownTags') && numel(im_info.UnknownTags.Value) == 1 %time tags from NSLS and APS large volume press beamlines -- if they exist
        im_times(i) = im_info.UnknownTags.Value;
    elseif isfield(im_info, 'UnknownTags') && numel(im_info.UnknownTags.Value) >= 1 %time tags from NSLS and APS large volume press beamlines -- if they exist
        im_times(i) = im_info.UnknownTags(1).Value;
    elseif strcmpi(type,'edf') == 1 %time stamps for ESRF 6ID image files
        dat = pmedf_read(im_name);
        loc = strfind(dat, 'time_of_day');
        im_times(i) = str2num((dat(loc+13:loc+31)));
    elseif strcmpi(type,'nxs') == 1 %time stamps for Diamond NXS image files
        im_times(i) = Image_Functions_NXSlist('times', im_name, 'camera', camera);
    else %if not read the fime the image was modified.
        im_times(i) = allimages(i).datenum * 24*60^2;
    end
end

% replace the time stamps from the save time with those in the replacement timestamp file (generally a log file).
if ts == 1
    
    if plot_on == 1
        im_t_old = im_times;
    end
    
    for i = 1:num_images
        loc = find(ismember(ImNames,  allimages(i).name));
        if isempty(loc)
            error_str = ['No time stamp has been found for ',allimages(i).name,' in the logfile.'];
            error(error_str);
        end
        im_times(i) = ImTimes(loc);
    end
    
elseif ts == 2 
    
    if numel(ImTimes) < num_images
        error('There are not enough time stamps in the replcement file for this process to succeed')
    end
    
    [~,lfil,e_str] = fileparts(ls(logfile));
    fprintf(' Replacing the image time stamps with those from %s.\n', [lfil,e_str])
    
    % it is necessary to assume that the first time stamps are very close
    % in time. If they are out by ~3600 seconds this is because the time
    % stamps have been adjusted by day light saving and the 3600 needs to
    % be added/subtracted from the times.    
    slack = 1200;
    if abs(mean(ImTimes) -  mean(im_times)) >= 3600 - slack %2 is arbitrary and needs to be checked
        if ImTimes(1) -  im_times(1) > 0
            im_times = im_times + 3600;
        else 
            im_times = im_times - 3600;
        end
    end
    
    %loop through the time stamps and match them to the log file times.
    wind = 70; %size of window to match time stamps over
    Locs = zeros(size(im_times));
    for x = 1:num_images
        
        if ~isempty(rep) && x <= rep %if we are just replacing the time stamps
            Locs(x) = x;
        else %match the differences. 
            
            w = x+wind-1;
            if w > numel(ImTimes)
                w = numel(ImTimes);
            end
            differences = repmat(im_times(x),1,w-x+1) - ImTimes(x:w);
            
            [~,l] = min(abs(differences));
            
            Locs(x) = x+l-1; %keep array of offsets 
        end
    end
  
    %make sure the array of time stamp locations does not skip a value and then have two the same.
    for x = 2:num_images-1
        if Locs(x) == Locs(x-1)+2 && Locs(x) == Locs(x+1)
            Locs(x) = Locs(x-1) + 1;
        end
    end
    if Locs(end-1) == Locs(end) && Locs(end)<length(ImTimes)
        Locs(end) = Locs(end)+1;
    end
    
    %make sure that the number of dropped frames always increases...
    done = 0;
    st = find(isfinite(Locs), 1, 'first');
    while done == 0
        Locs_old = Locs;
        for x = st:length(Locs)-1
            if ImTimes(Locs(x)) == ImTimes(end)% Locs(x) == length(ImTimes) %if we have got to the end of the log file times
                %do nothing
            elseif Locs(x) >= Locs(x+1)
                Locs(x) = Locs(x+1)-1;
%                 disp('conidiotn 1')
            end
            if x>=2 && Locs(x-1) < Locs(x) && Locs(x) > Locs(x+1)
                Locs(x) = Locs(x-1);
            end
        end
        if Locs_old(st:end) == Locs(st:end)
            done = 1;
        else
            if (0)plot_on == 1
                figure(4), hold off
                plot(1:length(Locs), Locs_old-(1:length(Locs)), '.-', 1:length(Locs), Locs-(1:length(Locs)), 'r.-')
                xlabel('Image number in sequence')
                ylabel('Dropped frames')
                title('Interating through dropped frames')
%                 pause(.1)
            end
        end
    end
           
    % and make new time stamp array.
    im_times_new = ImTimes(Locs);
        
    %make sure arrays are the same size
    if numel(im_times) == numel(im_times_new)
        if plot_on == 1
            im_t_old = im_times;
        end
        im_times = im_times_new;
    else
        error('The new time stamps are not the same size as those they are replaceing')
    end
    
end

%get temperatures of images (if required)
if add_temp == 1
    for i = 1:num_images
        try
            im_info = imfinfo([target_dir, allimages(i).name]);
        catch
            im_info = [];
        end
%         file_parts = FileTitleInformation(allimages(i).name);
        if isfield(im_info, 'UnknownTags') && numel(im_info.UnknownTags.Value) >= 1 %temperature from NSLS and APS large volume press beamlines -- if they exist
            im_temp(i) = str2num(im_info.UnknownTags(end).Value(13:end));
        elseif isfinite(file_parts(i,3)) %the image name has a temperature in it.
            im_temp(i) = file_parts{3};
        else %if not skip the image temperature.
            im_temp(i) = NaN;
        end
    end
end

%get power to images (if required)
if ~isempty(power_file) == 1
    p = importdata(power_file);
    pwt = datenum(p.textdata);
    pwr = p.data;
    for i = 1:num_images
        try
            im_info = imfinfo([target_dir, allimages(i).name]);
        catch
            im_info = [];
        end
        
        if isfield(im_info, 'FileModDate')
            im_time = im_info.FileModDate;
            im_time = datenum(im_time);
            
            %correct for US-UK time differences 
            im_time = im_time - 5/24; %(if needed).
        else
            error('Oh dear')
        end
        
        d = pwt-im_time;
        
        %plot(d)%, pause
        loc = find(min(abs(d))== abs(pwt-im_time));
        
        pow(i) = mean(pwr(loc));
            
    end
end


%sort the images by time. 
if numel(unique(im_times)) == numel(im_times)
    [im_times, order] = sort(im_times,direction);
else
% If the time stamps are not unique sort by image number
    for x = 1:length(allimages)
        aa(x,:) = FileTitleInformation(allimages(x).name, 'warn off');
    end
    [~, order] = sort(cell2mat(aa(:,5)),direction);
    im_times = im_times(order);
end
allimages = allimages(order);
if add_temp == 1
    im_temp = im_temp(order);
end

%make Lst filename
seps = find(target_dir == filesep);
if seps(end) == length(target_dir)
    file_name = target_dir(seps(end-1)+1:seps(end)-1);
else
    file_name = target_dir(seps(end)+1:end);
end

%catch the lst file name if it is something generic like 'image' and 
% replace it with something useful
if (~isempty(strfind(lower(file_name), 'image')) || ~isempty(strfind(lower(file_name), 'xrd'))) 
    
    expt = unique(file_parts(:,1));
    load = unique(cell2mat(file_parts(~isnan(cell2mat(file_parts(:,2))),2)));
    temp = unique(cell2mat(file_parts(~isnan(cell2mat(file_parts(:,3))),3)));
    per  = unique(cell2mat(file_parts(~isnan(cell2mat(file_parts(:,4))),4)));
    
    if numel(expt)==1 && numel(load)<=1 && numel(per)<=1 && numel(temp)<=1
        namGaps =  strfind(file_name,'_');
        if isempty(namGaps)
            namGaps = length(file_name);
        end
        gaps = strfind(allimages(1).name, '_');
        if isempty(gaps)
            gaps = length(allimages(1).name);
        end
        file_name = [allimages(1).name(1:gaps(end)-1),'_',file_name(1:namGaps(1))];
    end
end

%determine filename ending
if strcmpi(type, 'edf') == 1
    ending = 'edflst';
elseif strcmpi(type, 'nxs') == 1
    ending = 'nxslst';
else
    ending = 'lst';
end
file_name = [file_name, '.', ending];


%plot time stamps if requested
if plot_on == 1
    figure(1),
    clf
    hold on
    leg_str = [];
    if ts == 2 | ts==1
        plot(im_t_old, 1:length(im_t_old), 'go') %plot times images were saved at
        leg_str = [leg_str, {'Image Save times'}];
    end
    if ts ~= 0
        plot(ImTimes, 1:length(ImTimes), 'kx-') %plot times read from Log File.
        leg_str = [leg_str, {'Log File times'}];
    end
    plot(im_times, 1:length(im_times), '.-r') % plot time stamps to save
    leg_str = [leg_str, {'Time Stamps to save'}];
    xlabel('Time stamps')
    ylabel('Image number in sequence')
    title(['Times for ',file_name], 'Interpreter', 'None')
    legend(leg_str, 'Location', 'SouthEast')
    hold off
    
    if ts==2
        figure(2)
        clf
        hold on
        plot(1:length(Locs), Locs-(1:length(Locs)), '.-')
        xlabel('Image number in sequence')
        ylabel('Dropped frames')
        title('Dropped Frames')
        hold off
    end
    if ts~=0        
        figure(3)
        clf
        hold on
        plot(1:length(im_times), im_times-im_t_old, 'b.-')
        xlabel('Image number in sequence')
        ylabel('Time difference between save time and log file time')
        title('Delta T')
        hold off
    end
    
    pause(1)
end


%write list of file names and timestamps to the lst file.
outfile = fopen(file_name, 'w');
for x = 1 : num_images
    str1 = [target_dir,allimages(x).name];
 
%     fprintf(outfile, '%s, %20.5f,\n', str1, im_times(x));
    if add_time == 0
        fprintf(outfile, '%s,\n', str1);
    elseif add_temp == 0
        fprintf(outfile, '%s, %20.5f,\n', str1, im_times(x));
    elseif ~isempty(power_file) && add_temp == 0
        fprintf(outfile, '%s, %20.5f, , %12.3f\n', str1, im_times(x), pow(x));
    elseif ~isempty(power_file) && add_temp == 1
        fprintf(outfile, '%s, %20.5f, %12.3f, %12.3f\n', str1, im_times(x), im_temp(x), pow(x));
    else
        fprintf(outfile, '%s, %20.5f, %12.3f,\n', str1, im_times(x), im_temp(x));
    end
        
end

fclose(outfile);

%% versions/changes
% v 1.4.2 -- 18th Feb 2019
%   - fixed error in inputting camera type     
% v 1.4.1 -- 26th November 2018
%   - added *.db as a disallowed file ending for the listings.
% v 1.4 -- 15th October 2018
%   - Added ability to make list file from nxs files (from Diamond).
%   - list file conatining nxs files ends with .nxslst
% v 1.3 -- 3rd July 2018
%   - Added options so that can include power in the list files..
% v 1.2 -- 21st August 2017
%   - Edited so that can write a list file without the time stamps.
% v 1.1 -- 27 June 2017
%	- Rewrote replace function so that it works better. It now works on
%	much poorer time stamp arrays.
%   - Added plot function to plot time stamps before ending the function.
% v 1.0 -- 23 June 2017
%	- implemented log file and replacemnt time stamps. + fixed some bugs therein.
% 9/7/17
% - started implementing log file timestamps.
% 8/4/17
% - added option to read temperatures from tags in file and list with times
%   This is for the ImageAnalysisLoop scripts. 
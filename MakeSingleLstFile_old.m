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
%
%    'image format': recognsed image formats are 'tif', 'tiff', 'edf' and
%                    'bmp'. If none is specified 'tif' is the default.
%    'order':        'ascend' or 'descend' to determine the order the images are
%                    listed in. 'ascend' is the default.
%    'temperature':  Option to read tempeartures from image files and list in the lst file
%    'logfile':      Option to read time stamps from a logfile.
%    'replace':      Option to read time stamps from another file, this just replaces the first 
%                        Nkeep time stamps by postion in the list and then minimises the differences 
%                        for the rest.
%                   fname can take the format asdf->fgjk/*.log -- in which
%                   case the script replaces asdf with fghk in the target
%                   directory string and appends everything after the file
%                   separator.
%
%    For both logfile and replace options a relative location is assumed for the log file 
%       if the string starts with './', '.\', '../' or '..\'.
%
%   See Also MakeManyTimes, MakeManyLstFiles

%	Version 1.0
%   Simon Hunt 2016-2017
% This script used to be called makemanytimes_tiff.

function  MakeSingleLstFile(target_dir, varargin)

%process varargin and defaults
plot_on = 0;
direction = 'ascend';
type = 'tif';
add_temp = 0;
logfile = [];
ts = 0; %time stamp replacement settings: 0 = void; 1 = logfile [match file names with log file]; 2 = replace [match save time with closest time in logfile]
rep = [];

iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case {'descend', 'ascend'}
            direction = varargin{iarg};
            iarg = iarg + 1;
        case {'tiff', 'tif', 'bmp', 'edf'}
            type = varargin{iarg};
            iarg = iarg + 1;
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
            if isnumeric(varargin{iarg})
                rep = varargin{iarg};
                iarg = iarg + 1;                
            end
%             error('Logfile input is not implemented.')
        case {'temperature'}
            add_temp = 1;
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
    endings(any_dir) = []; %remove 'dir' as an option
    endings(any_lst) = []; %remove '.lst' as an option
    endings(any_chi) = []; %remove '.lst' as an option
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
num_images = length(allimages);

fprintf(' Found %i *.%s files to add to .lst file.\n', num_images, strrep(type, '.', ''))

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
    try
        im_info = imfinfo([target_dir, allimages(i).name]);
    catch
        im_info = [];
    end
    if ts == 1 %get the times of the images from the logfile.        
        loc = find(ismember(ImNames,  allimages(i).name));
        if isempty(loc)
            error_str = ['No time stamp has been found for ',allimages(i).name,' in the logfile.'];
            error(error_str);
        end
        im_times(i) = ImTimes(loc);       

    elseif isfield(im_info, 'UnknownTags') && numel(im_info.UnknownTags.Value) == 1 %time tags from NSLS and APS large volume press beamlines -- if they exist
        im_times(i) = im_info.UnknownTags.Value;
    elseif isfield(im_info, 'UnknownTags') && numel(im_info.UnknownTags.Value) >= 1 %time tags from NSLS and APS large volume press beamlines -- if they exist
        im_times(i) = im_info.UnknownTags(1).Value;
    elseif strcmpi(type,'edf') == 1 %time stamps for ESRF 6ID image files
        dat = pmedf_read([target_dir,allimages(i).name]);
        loc = strfind(dat, 'time_of_day');
        im_times(i) = str2num((dat(loc+13:loc+31)));
    else %if not read the fime the image was modified.
        im_times(i) = allimages(i).datenum * 24*60^2;
    end
end

% replace the time stamps from the save time with those in the replacement timestamp file (generally a log file).
if ts == 2 
    
    if numel(ImTimes) < num_images
        error('There are not enough time stamps in the replcement file for this process to succeed')
    end
    
    [~,lfil,e_str] = fileparts(ls(logfile));
    fprintf(' Replacing the image time stamps with those from %s.\n', [lfil,e_str])
    
    % it is necessary to assume that the first time stamps are very close
    % in time. If they are out by ~3600 seconds this is because the time
    % stamps have been adjusted by day light saving and the 3600 needs to
    % be added/subtracted from the times.    
    slack = 600;
    if abs(ImTimes(1) -  im_times(1)) >= 3600 - slack %2 is arbitrary and needs to be checked
        if ImTimes(1) -  im_times(1) > 0
            im_times = im_times + 3600;
        else 
            im_times = im_times - 3600;
        end
    end
    
    wind = 70; %size of window to match time stamps over
    for x = 1:num_images
        
        if ~isempty(rep) && x <= rep %if we are just replacing the time stamps
            im_times_new(x) = ImTimes(x);
            Locs(x) = NaN;
        else %match the differences. 
            
            w = x+wind-1;
            if w > numel(ImTimes)
                w = numel(ImTimes);
            end
            differences = repmat(im_times(x),1,w-x+1) - ImTimes(x:w);
            
            [a,l] = min(abs(differences));
            
            Locs(x) = l-1; %keep array of offsets 
        end
    end
    
%     % subtle check for dropped frames at high frame rate -- check time
%     % diference averages
% %     ave_over = 3;
%     done = 0;
%     while done == 0
%         t_diff = im_times-ImTimes(1:length(im_times));
%         t_diff = t_diff(:);
%         ave = mean([t_diff(1:end-2), t_diff(2:end-1), t_diff(3:end)],2);
%         dif_ave = ave(rep:end-1) - ave(rep+1:end);
%         bigbest = find(abs(dif_ave) == max(abs(dif_ave)));
%     end
        

    %make sure the differences do not swtich between being positive and
    %negatvie
    d = 0;
    it_pos = rep:length(im_times);
    while d == 0 
        diffs = (ImTimes(rep:length(im_times))) - im_times(it_pos);
        
        small = find(diffs<0);
        sm = find(small~=1, 1, 'first')
       
        if isempty(sm)
            d = 1;
        else
            it_pos(small(sm):end) = it_pos(small(sm):end)+1;
            
            plot(it_pos,'.')
            keyboard
        end

    end
    
    % make sure the array of offests does not bounce between values (i.e. 1,0,1,0,1,0...) 
    for x = 2:num_images-1
        if Locs(x) == Locs(x-1) + 1 && Locs(x) == Locs(x+1) + 1
            Locs(x) = Locs(x) - 1;
        end
    end
    
    %make sure that the number of dropped frames always increases...
    done = 0;
    st = find(isfinite(Locs), 1, 'first');
    while done == 0
        Locs_old = Locs;
        for x = st:length(Locs)-1
            if Locs(x) > Locs(x+1)
                Locs(x) = Locs(x+1);
            end
            if x>=2 && Locs(x-1) < Locs(x) && Locs(x) > Locs(x+1)
                Locs(x) = Locs(x-1);
            end
        end
        if Locs_old(st:end) == Locs(st:end)
            done = 1;
        else
            if plot_on == 1
                plot(1:length(Locs), Locs_old, '.-', 1:length(Locs), Locs, 'r.-')
                xlabel('Image number in sequence')
                ylabel('Dropped frames')
                pause(.1)
            end
        end
    end
           
    % and make new time stamp array.
    for x = 1:num_images
        if isnan(Locs(x))
            im_times_new(x) = ImTimes(x);
        else
            im_times_new(x) = ImTimes(x+Locs(x));
        end
    end
    
    %plot changes in required
    if plot_on == 1
        figure, 
        hold on
        plot(ImTimes-im_times(1), 1:length(ImTimes), 'k.-')
        plot(im_times-im_times(1), 1:length(im_times), '.-r')
        plot(im_times_new-im_times(1), 1:length(im_times_new), 'g.')
        legend('Log file times', 'Image save times', 'Replacement image times', 'Location', 'SouthEast')
        xlabel('Time stamps')
        ylabel('Image number in sequence')
        
        
        figure 
        hold on
        plot(1:length(Locs), Locs, '.-')
        xlabel('Image number in sequence')
        ylabel('Dropped frames')
        
        keyboard
    end
    
    %make sure arrays are the same size
    if numel(im_times) == numel(im_times_new)
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
% 
% manual = 1;
% 
% if manual == 1
%     d = 0;
%     while d == 0
%         
%         rep = input('Additional dropped images? (one at a time)')
%     
        
    
% %sort the images by time. 
% if numel(unique(im_times)) == numel(im_times)
%     [im_times, order] = sort(im_times,direction);
% else
% % If the time stamps are not unique sort by image number
%     for x = 1:length(allimages)
%         aa(x,:) = FileTitleInformation(allimages(x).name, 'warn off');
%     end
%     [~, order] = sort(cell2mat(aa(:,5)),direction);
%     im_times = im_times(order);
% end
% allimages = allimages(order);
% if add_temp == 1
%     im_temp = im_temp(order);
% end

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
    
    if numel(expt)==1 & numel(load)<=1 & numel(per)<=1 & numel(temp)<=1
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


%write list of file names and timestamps to the lst file.
outfile = fopen([file_name,'.lst'], 'w');
for x = 1 : num_images
    str1 = [target_dir,allimages(x).name];
 
%     fprintf(outfile, '%s, %20.5f,\n', str1, im_times(x));
    if add_temp == 0
        fprintf(outfile, '%s, %20.5f,\n', str1, im_times(x));
    else
        fprintf(outfile, '%s, %20.5f, %12.3f,\n', str1, im_times(x), im_temp(x));
    end
        
end

fclose(outfile);

%% versions/changes
% v 1.0 -- 23 June 2017
%	- implemented log file and replacemnt time stamps. + fixed some bugs therein.
% 9/7/17
% - started implementing log file timestamps.
% 8/4/17
% - added option to read temperatures from tags in file and list with times
%   This is for the ImageAnalysisLoop scripts. 
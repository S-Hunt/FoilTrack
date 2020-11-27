% MakeManyLstFiles(directory,...).
%
% Makes a *.lst file for the images in the target directory.
% In the list file is the list of images and their associated time stamps. 
% It is assumed that the data images are the most common in the target directroy. 
%
% Syntax:
%       MakeManyLstFiles('target directory')
%       MakeManyLstFiles('target directory', 'subdirs', 0/1/2)
%       MakeManyLstFiles('target directory', 'image format')
%       MakeManyLstFiles('target directory', 'order')
%       MakeManyLstFiles('target directory', 'temperature')
%       MakeManyLstFiles('target directory', 'logfile', 'fname')
%       MakeManyLstFiles('target directory', 'replace', fname, [Nkeep])
%       MakeManyLstFiles('target directory', 'plot')
%
%    'subdirs':      not required but if present (default is 2): 
%                      0 - only makes list file for target directory
%                      1 - only makes list file for subdirectories of the target
%                      2 - makes list files for directory and its subdirectories recursively
%                   'string' - makes list files for only the directories
%                              with the strain in their directory path.
%    'image format': recognsed image formats are 'tif', 'tiff', 'edf', 'netcdf', 'nc' and
%                    'bmp'. If none is specified MakeSingleListFile determines the format.
%    'order':        'ascend' or 'descend' to determine the order the images are
%                    listed in. 'ascend' is the default.
%    'temperature':  Option to read tempeartures from image files and list in the lst file
%    'logfile':      Option to read time stamps from a logfile.
%    'replace':      Option to read time stamps from another file, this just replaces the first 
%                        Nkeep time stamps by postion in the list and then minimises the differences 
%                        for the rest.
%    'plot':          Plot the time stamps vs. frame number for each lst file.
%
%    For both logfile and replace options a relative location is assumed for the log file 
%       if the string starts with './', '.\', '../' or '..\'.
%
%   See Also MakeManyTimes, MakeSingleLstFile

%   Simon Hunt 2007 - 2018
%   Version 1.3

function  MakeManyLstFiles(target_dir, varargin)

%process varargin and defaults
direction = 'ascend';
image_type = [];
subdirs = 2;
add_temp = 0;
power_file = [];
lf = 0;
logfile = [];
rep = [];
plot_on = 0;

iarg = 1;
while iarg <= (length(varargin))
    switch lower(varargin{iarg})
        case {'descend', 'ascend'}
            direction = varargin{iarg};
            iarg = iarg + 1;
        case {'tiff', 'tif', 'bmp', 'netcdf', 'nc', 'nxs', 'edf'}
            image_type = varargin{iarg};
            iarg = iarg + 1;
        case {'subdirs'}
            subdirs = varargin{iarg+1};
            iarg = iarg + 2;
        case {'logfile'}
            lf = 1;
            logfile = varargin{iarg+1};
            iarg = iarg + 2;
        case {'replace'}
            lf = 2;
            logfile = varargin{iarg+1};
            iarg = iarg + 2;
            if isnumeric(varargin{iarg})
                rep = varargin{iarg};
                iarg = iarg + 1;                
            end
        case {'temperature'}
            add_temp = 1;
            iarg = iarg + 1;
        case {'power'}
            power_file = varargin{iarg+1};
            iarg = iarg + 2;
        case 'plot'
            plot_on = 1;
            iarg = iarg + 1;
        otherwise
            error(['Unknown option: ' varargin{iarg}]);
    end
end

%check trget directory exists and is correctly formatted.
if exist(target_dir, 'dir') ~= 7
    error('The target folder does not exist');
end

if target_dir(end) ~= '/' && target_dir(end) ~= '\'
    target_dir = strcat(target_dir,filesep);
end

%arguemnts to pass to MakeSingleLstFile
args_to_pass = {direction};
if ~isempty(image_type) %pass image type if set
   args_to_pass = [args_to_pass, {image_type}]; 
end
if lf == 1 % pass logfile or replce settings if called
   args_to_pass = [args_to_pass, {'logfile'}, {logfile}]; 
elseif lf == 2 && isempty(rep)
   args_to_pass = [args_to_pass, {'replace'}, {logfile}]; 
elseif lf == 2 && ~isempty(rep)
   args_to_pass = [args_to_pass, {'replace'}, {logfile}, {rep}]; 
end
if add_temp == 1 %pass temperature call. 
   args_to_pass = [args_to_pass, {'temperature'}];
end
if ~isempty(power_file) %pass power call. 
   args_to_pass = [args_to_pass, {'power', power_file}];
end
if plot_on == 1
    args_to_pass = [args_to_pass, {'plot'}];
end
    

%% make lst file for subdirectories

if ischar(subdirs) | ~isequal(subdirs, 0)
    
    %get listing of all subdirs.
    if ispc
        listing = regexp(genpath(target_dir),'[^;]*','match')'
    else
        listing = regexp(genpath(target_dir),'[^:]*','match')'
    end
        
    %cut list to only directories of interest 
    if ischar(subdirs)
        loc = strfind(listing, subdirs);
        
        remove = [];
        for x = 1:length(listing)
            if isempty(loc{x})
                remove = [remove, x];
            end
        end
        listing(remove) = [];
    end
    
    for x = 1:length(listing)
        listing_target = listing{x};
        
        fprintf('Making .lst file for %s\n',listing_target);
        MakeSingleLstFile(listing_target,args_to_pass{:}); %make lst file in sub directory
    end
end

%% make lst files for current directory.
% only required if listing for current directory only
if isequal(subdirs, 0)
    fprintf('Making .lst file for %s\n',target_dir);
    MakeSingleLstFile(target_dir,args_to_pass{:}); %make lst file in current directory
    
end


%% Versions
% - 1.3.1 - 15th October 2018
%   - added option to include 'nxs' and 'edf' images as options.
% - 1.3 - 3rd July 2018
%   - added option to include power in the list files..
% - 1.2.2 - 27rd June 2017
%   - added function to pass plots to MakeSingleLstFile.
% - 1.2.1 - 23rd June 2017
%	- some bug fixes so that script runs properly on linux/r2015a
% - 1.2 - 9 to 13th June 2017
%   - Added option to restrict listing to directories with sting in file
%   path.
%   - Added options to read time stamps from a log file and/or replace time
%   stamps with those from log file (for anelasticity in DAC expeirments
%   with crap time stamps).
% - 1.1.2 - 18 May 2017
%   - Changes so that the script works in linux. Specifically so that the
%   listing of all the subdirectories works.
% - 1.1.1 - 8 april 2017
%   - Added option for temperatures in the list files.
% - 1.1 - 16th Jan 2017
%   - Changed the directory search function so that it recursively searches all subfolders.
% - 1.0 
%   - all changes prior to end of 2016
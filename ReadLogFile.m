% ReadLogFile(fname).
%
% Reads the log files to extract time stamps. Called by MakeSingleLstFile. 
%
% Syntax:
%       [times, names] = ReadLogFile(fname)
%
% The script assumes that the logfile is of the format:
%      endTime               =========== File ============
%   yyyy-mm-ddThh:mm:ss.sss  /dir/ImageName1.ext
%   yyyy-mm-ddThh:mm:ss.sss  /dir/ImageName2.ext
%   yyyy-mm-ddThh:mm:ss.sss  /dir/ImageName3.ext
%   yyyy-mm-ddThh:mm:ss.sss  /dir/ImageName4.ext
%   .....
%   0.4000  0.4000  2016-11-03T01:39:52.878  /dir/ImageNameN.ext
%   --------- End ----------------------------
%
%  This is the log file format for the xrd images from APS, 16IDB in Nov
%  2016. The experiments with the data was BCC1/HCP1.
%
%   See Also MakeManyTimes, MakeManyLstFiles, MakeSingleLstFile

%   Simon Hunt 2017

function [TimeStamp names] = ReadLogFile(fname)


%check there is only one file denoted by 'fname'
f = dir(fname);

if numel(f) > 1
    error('there is more than one file denoted by the logfile string')
elseif isempty(f)
    error('the logfile cannot be found')
end

%this is needed to allow for the fname to refer to a wildcard which fopen
%does not like.
[d,~,~]=fileparts(fname);

fid = fopen([d,filesep,f.name]);
%importdata(fname);


done = 0;
line = 0;

TimeStamp = [];
names = [];

while done == 0
    
    tline = fgetl(fid);
    
    if tline == -1
        done = 1;
    else
        
        str = textscan(tline, '%s');
        line = line+1;
        
        if line == 1 %ignore the header line.
           
        elseif strcmpi(str{:}{1}, '###') == 1 %ignore lines starting with ### 
            
        elseif numel(str{:}) == 2 || numel(str{:}) == 4 % list the file names and times.
            if numel(str{:}) == 2 
                loc = 1;
            elseif numel(str{:}) == 4
                loc = 3;
            end           
            
            %get the time stamp
            TimeStr = str{1}{loc};
            TimeStr = strrep(TimeStr, 'T', ' '); %This is to replce the T in the log files from APS 16IDBxrd (BCC1 anelasticity of Fe experiments).
            
            %parse time string
            TimeStamp(end+1) = datenum(TimeStr) * 24*60^2; %convert date into seconds
            

            %get the file name
            %keep the file name but ditch the directory.
            [~,nam,e] = fileparts(str{1}{loc+1});
            names{end+1} = [nam,e];            
            
        end
        
    end
    
    
end
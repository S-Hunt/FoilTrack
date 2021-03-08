%% wildcard2filelist.m
%% Makes a list of all files given by a shell wildcard, e.g. 'Si_*.edf' or
%% '*.txt'. Files are obtained by the "ls" command, and their names are
%% returned in the resulting files{1..n} cell array.
%%
%% Syntax:
%%	files = wildcard2filelist( wild )
%%
%% Example 1:
%%	files = wildcard2filelist('*.txt');
%%	fprintf('Found %i such files, first is %s\n', length(files), files{1});
%%
%% Limitations: don't use it for something like '-R /usr/local/lib' -- there,
%% another version of this program would probably be more useful -- such one
%% using the 'find' program.
%%
%% Author: Petr Mikulik
%% Version: April 2009
%% History:
%%     1.4.2009: "\n" changed to sprintf('\n') for Matlab compatibility
%%    16.6.2004: Original version.

function files = wildcard2filelist( wildcard )

if nargin~=1
    error('Syntax: files = wildcard2filelist( wild )');
end

% Get list of files by the 'ls' command:
[res,ls] = unix(['ls -1 ', wildcard, ' 2>/dev/null']);

n = sum(ls==sprintf('\n'));
if n==0
    % fprintf(['WARNING: No file named as "', wildcard, '" found.\n']);
    files={};
    return
end

k=0;
while 1
    [name, ls]=strtok(ls,sprintf('\n'));
    if length(name)==0 break; end
    k=k+1;
    files{k}=name;
end

% eof wildcard2filelist.m

%% matbinning.m
%%
%% Syntax:
%%	b = matbinning ( a, binning )
%% where a is vector or matrix and binning is vector of binnings rows and cols.
%%
%% Example:
%%	matbinning(1:10, [1,2])
%%	matbinning(1:10, [1,3])
%%	matbinning(1./hilb(10), [2,5])
%%
%% Autor: Petr Mikulik
%% Version: August 2010
%% Previous version: September 2004

function a = matbinning ( a, binning )

if nargin~=2
    error('Usage: b = matbinning(a, [binrows, bincolumns])')
end
if length(size(a))~=2
    error('a must be matrix');
end
binning=binning(:)';
if length(binning)~=2
    error('binning must have same size as a (e.g. 2x2)');
end


if any(binning > 1) % do binning of a matrix
    [nr, nc]=size(a);
    if nr>1 nr=binning(1)*floor(nr/binning(1)); end
    if nc>1 nc=binning(2)*floor(nc/binning(2)); end
    a=a(1:nr, 1:nc); % cut extra elements
    if binning(1)>1
	a=reshape(a, binning(1), round(nr*nc/binning(1)));
	a=reshape(sum(a), nr/binning(1), nc);
	nr=floor(nr/binning(1));
    end
    if binning(2)>1
	a=reshape(a', binning(2), round(nr*nc/binning(2)));
	a=reshape(sum(a), nc/binning(2), nr)';
    end
end

% eof matbinning.m

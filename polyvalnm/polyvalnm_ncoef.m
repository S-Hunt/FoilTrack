% polyvalnm_ncoef(n,m).
% Calculates the expected number of coeficients for a n,m polynomial
% surface.
%
% The order of the surface is max(n, m).
%
% See also: Polyvalnm, Polyvalnm_coef2v, Polyvalnm_coef2mat

% Simon Hunt 21 May 2014, 19th Dec 2017.

function out = polyvalnm_ncoef(m,varargin)

%parse input
if numel(m) == 2 && isempty(varargin)
    n = m(2); %these are delibarately not alphabetical. 
    m = m(1);
elseif numel(m) == 1 && nargin == 2
    m = m;
    n = varargin{1};
else
    error('The number of inputs is wrong')
end   
    
%process
M = m+1;
N = n+1;

if n==0 || m==0
    out = max(m,n) + 1;
else
    out = M*N - nchoosek(min(M,N),2);
end
% Reshapes coefficients for polnomial surface of order n,m into matrix
% Where n and m are the orders of the x and y polnomials 
% void numbers are set as NaN

function out = polyvalnm_coef2mat(p,n,m,varargin)

%defines dimensions of (re)shaped array
M = m+1;
N = n+1;

%determines what the void values in the matrix should be.
if length(varargin) >= 1
    if strcmpi(varargin{1},'zeros')
        type = 'zero';
    else
        type = 'NaN';
    end
else
    type = 'NaN';
end

%Checks length of input vector
q = p(:);
q(isnan(q)==1) = [];
if (length(q) ~= polyvalnm_ncoef(n,m) && isnan(p(1))) || (length(q) ~= N*M && ~isnan(p(1)) && ~(size(p,1)== 1 || size(p,2)==1))
    error('Wrong number of coefficients')
end

if ~isequal(size(p),[N M]) 
    %turn linear array of coefficients into matrix n+1 * m+1 matrix with
    %NaN for unused coeficients
    [i j] = ind2sub([N,M],1:N*M);
    k = N*M:-1:1;
    for x = 1:length(k)
        if i(x)+j(x)>max(N,M)+1
            k(x) = NaN;
        end
    end
    k(isnan(k) == 1) = [];
    
    %create blank matrix
    if strcmpi(type, 'NaN')
        l = NaN(N,M);
    elseif strcmpi(type, 'zero')
        l = zeros(N,M);
    else
        error('Unknown type')
    end
    
    %populate values of p into matrix.
    l(k) = flipud(p(:));
    out = l;
    
elseif (strcmpi(type, 'zero') == 1 && isnan(p(1,1)))
    %change NaN values in p to zeros.
    p(isnan(p)) = 0;
    out = p;
    
elseif (strcmpi(type, 'NaN') == 1 && p(1,1) == 0)
    %change void values from 0 to NaN. 
    % This cannot be done the same as turning NaN to 0 incase any of the
    % parameters of the array are 0.
    [i j] = ind2sub([N,M],1:N*M);
    k = N*M:-1:1;
    for x = 1:length(k)
        if N-i(x)+1+M-j(x)+1>max(N,M)+1
            p(x) = NaN;
        end
    end   
    
    out = p;
    
else
    %if p is already the size of [N M] then just return it.
    out = p;
    
end
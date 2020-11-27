% Reshapes coefficients for polnomial sirface of order n,m
% Where n and m are the orders of the x and y polnomials 

function out = polyvalnm_reshape(p,n,m)

%defines dimensions of (re)shaped array
M = m+1;
N = n+1;

%Checks length of input vector
q = p(:);
q(isnan(q)==1) = [];
if length(q) ~= polyvalnm_ncoef(n,m)
    error('Wrong number of coefficients')    
end
     
     
if length(p) == numel(p) %FIXME: this does not work for 'normal' polynomials where n or m = 0
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
    l = NaN(N,M);
    l(k) = flipud(p);
    out = l;
    
else  
    %turn matrix p into vector with no NaN values in it.
    p = p(:);
    p(isnan(p)==1) = [];
    out = p;
end
% Reshapes coefficients for polnomial surface of order n,m to a vector.

function out = polyvalnm_coef2v(p)


if ~isnan(p(1)) && (size(p,1)~=1 || size(p,2)~=1)
    
    p = polyvalnm_coef2mat(p,size(p,1)-1,size(p,2)-1,'NaN')
    
end
    
    %turn matrix p into vector with no NaN values in it.
    p = p(:);
    p(isnan(p)==1) = [];
    out = p;

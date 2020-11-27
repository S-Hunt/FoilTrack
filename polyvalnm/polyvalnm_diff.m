% polyvalnm_diff(p,dim,n)
% dim - dimension, n - nth root, 
% Calculates the differential of the 2D polnomial surface fed in. 
% of order max(n,m)


function out = polyvalnm_diff(p,dim,n,varargin)

if dim ~= 1 && dim ~= 2
    error('Unknown dimension to differentiate in.')
end

% reorient the coefficient matrix.
if dim == 2
    p = p';
end

for x = 1:n
    order = size(p);
    [~, Y]=meshgrid(0:order(2)-1, order(1)-1:-1:0);
    p = p.*Y;
    
%     p(1,:) = [];
    p(end,:) = []; %the above line was wrong. I think that this is correct - 20 Dec '17
    if isnan(p(:,end)) == 1  % these lines are probably unnecessary.
        p(:,end) = [];
    end
    if isnan(p(:,1)) == 1
        p(:,1) = [];
    end
end

if dim == 2
    p = p';
else
    p = p;
end
out = p;


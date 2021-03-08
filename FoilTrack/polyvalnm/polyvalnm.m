%calculate values for 2D polnomial surface of order max(n,m)


function out = polyvalnm(p,x,y,varargin)

size_in = size(x);
if size(x,1) ~= size(y,1) | size(x,2) ~= size(y,2)
    error('The x and y vectors must be the same size')
else
    x = x(:);
    y = y(:);
end

%switch orientations if p is deeper than it is wide.
if size(p,1) > size(p,2)
   p = p';
   x_in = x;
   x = y;
   y = x_in;
   clear x_in
end


N = size(p,1);
M = size(p,2);

if length(varargin) >= 2
    n = varargin{1};
    m = varargin{2};
else
    n = N-1;
    m = M-1;
end

if numel(p) ~= N*M && (size(p,1)~=1 || size(p,1)~=1)
    error('Wrong number of coeficients') 
end

p = polyvalnm_coef2mat(p,n,m,'NaN');


% [X,Y] = meshgrid(M-1:-1:0, N-1:-1:0);
% X(isnan(p(:))) = NaN;
% Y(isnan(p(:))) = NaN;

% p = polyvalnm_coef2v(p);
% X = polyvalnm_coef2v(X);
% Y = polyvalnm_coef2v(Y);

a = 0;
for i = 1:N
    
    coeffs = p(i,:);
    coeffs(isnan(coeffs)==1) = [];

%     a(i,:) = polyval(coeffs,x) .* y.^(N-i);
% polyval(coeffs,x) .* y.^(N-i)
    a = a + polyval(coeffs,x) .* y.^(N-i);
    
end

out = reshape(a,size_in);


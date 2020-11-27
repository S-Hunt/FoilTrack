% polyvalnm_solve(p, equal_to, x_vals, [y_range])
% p - polynomial, equal_to - value to solve for, x_vals -- for these x values, 
% Solves the polynomial for p(x_vals) = equal_to. 
% Calculates the polynomial and then fits a new polynomial to the values.


function out = polyvalnm_solve(poly, equal_to, x_vals, varargin)

% y -- displacement
% x -- time.
% solving for min x at each y.

if nargin >=4
    range = varargin{1};
else
    range = [-10, 10];
end

%FIX ME: assumes the polynomial coefficients are a matix
%FIX ME: Turns equal_to into a string before solving.

%p = polyvalnm_coef2mat(p);

[n m] = size(poly);
n = n-1;
m = m-1;

if n == 0 || m == 0
    
   %
   disp('Need to solve for single array not a matrix')
   return
end
  
p_keep = poly(:);%polyvalnm_coef2v(p);
    
poly_str = polyvalnm_symbolic(polyvalnm_coef2v(poly),n,m);


syms f(x,y,p)

eqn = [poly_str, ' == ', num2str(equal_to)];

eqn_sub1 = eqn;

for j = 1:numel(p_keep);
    eqn_sub1 = strrep(eqn_sub1, ['p(',num2str(j),')'], num2str(p_keep(j)));
end

a = NaN(1,length(x_vals));

for j = 1:length(x_vals)
    
     
    %print running tally of where we are. 
    str = sprintf('  Solving point %i of %i', j, length(x_vals));
    del = sprintf(repmat('\b',1,length(str)));
    fprintf(str)
    
%      eqn_sub = subs(eqn_sub1, y, y_vals(j),0);
     eqn_sub = subs(eqn_sub1, x, x_vals(j));%,0);
   
     if eqn_sub == 1 || eqn_sub == 0
          nums = solve(eqn_sub1);
     else
         vals = solve(eqn_sub, y);
         
         nums = eval(vals);
     end
         
     %discard complex roots.
     for k = 1:numel(nums)
         re(k) = isreal(nums(k));
     end
     nums(re==0) = [];
     
     %cut roots to be within range.
     lo = nums<min(range);
     nums(lo) = [];
     hi = nums>max(range);
     nums(hi) = [];
     
     %make sure only have one copy of each root.
     nums = unique(nums);
     
     if numel(nums) == 1
         a(j) = nums;
     elseif isempty(nums)
         a(j) = NaN;
     else
         a(j) = NaN;
     end
     
     fprintf(del)
end  
     
out = a;

% 
% return
% 
% 
% sol = solve(eqn,y);
% 
% f(x,y,p) = sol
% 
% solve(f,y, p, p_keep', x , 0)
% 
% syms p x y
% eqn = [poly_str, ' = ', num2str(equal_to)];
% sol = solve(eqn,y);
% 
% subs(solx, p, p_keep)
% 
% % out = solve(sol, num2str(x_vals), )
% 
% for j = 1:length(x_vals)
%     
%     y_min(j) = solve([sol, ' = ' num2str(x_vals(j))] );
%     
% end
%     
% 
% return
% 
% if dim ~= 1 && dim ~= 2
%     error('Unknown dimension to differentiate in.')
% end
% 
% % reorient the coefficient matrix.
% if dim == 2
%     p = p';
% end
% 
% for x = 1:n
%     order = size(p);
%     [~, Y]=meshgrid(0:order(2)-1, 0:order(1)-1);
%     p = p.*Y;
%     
%     p(1,:) = [];
%     if  isnan(p(:,end)) == 1
%         p(:,end) = [];
%     end
% end
% 
% if dim == 2
%     p = p';
% else
%     p = p;
% end
% out = p;


    
    
    
    




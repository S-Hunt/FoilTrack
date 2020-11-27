function Polyvalnm_testcases2
% POLYVALNM_TESTCASES. Runs set of test cases to show functions in the polyvalnm suite 
% are working correctly



%FIX ME.
%Need more test cases to:
% - confirm the coeficient making code gives correct outputs
% - the polyvalnm works with both 'zero' and 'NaN' options.

fprintf('\nPOLYVALNM_TESTCASES \n')

%% Test polyvalnm_ncoef
fprintf('Polyvalnm_ncoef, case 1: Number of coefficients for single row vector; numel(x) > 1, numel(y) = 1\n')

order = round(rand(1)*25);
 fprintf('    Order of polynomial: [0, %i]; expected number of coefficients = %i\n', order, order+1)

 fprintf('    Call: Polyvalnm_ncoef(0, %i)\n', order)
v = polyvalnm_ncoef(0,order);
 fprintf('    Number of coefficients from Polyvalnm_ncoef: %i\n', v)

pass = order+1 == v;

if pass == 1
 fprintf(' PASSED. Number of coefficients is correct\n\n')
elseif pass ~= 1
    keyboard
end

%case 2
fprintf('Polyvalnm_ncoef, case 2: Number of coefficients for single column vector; numel(x) = 1, numel(y) > 1\n')

order = round(rand(1)*25);
 fprintf('    Order of polynomial: [%i, 0]; expected number of coefficients = %i\n', order, order+1)

 fprintf('    Call: Polyvalnm_ncoef(%i, 0)\n', order)
v = polyvalnm_ncoef(order,0);
 fprintf('    Number of coefficients from Polyvalnm_ncoef: %i\n', v)

pass = order+1 == v;

if pass == 1
 fprintf(' PASSED. Number of coefficients is correct\n\n')
elseif pass ~= 1
    keyboard
end

%case 3
fprintf('Polyvalnm_ncoef, case 3: Number of coefficients for non singular array; numel(x) > 1, numel(y) > 1\n')

order = round(rand(1,2)*25)+1;
remove = min(order(1), order(2));
expectd_size = (order(1)+1)*(order(2)+1) - nchoosek(min(order)+1,2);
% expectd_size = 0;
% for x = 1:order(2)
%     expectd_size = expectd_size + order(1) - x+1;
% end
 fprintf('    Order of polynomial: [%i, %i]; expected number of coefficients = %i\n', order(1), order(2), expectd_size)

 fprintf('    Call: Polyvalnm_ncoef(%i, %i)\n', order(1), order(2))
v = polyvalnm_ncoef(order(1),order(2));
 fprintf('    Number of coefficients from Polyvalnm_ncoef: %i\n', v)

pass = expectd_size == v;
 
if pass == 1
 fprintf(' PASSED. Number of coefficients is correct\n\n')
elseif pass ~= 1
 fprintf(' FAILED. Number of coefficients is not the same\n\n')
    keyboard
end


fprintf('Calculate polynomials using Polyvalnm and compare to known functions.\n\n')

%% Case 1: coefficients are single row vector 
disp('Case 1a: Polynomial coefficients are single row vector; numel(x) > 1, numel(y) = 1')

l = round(rand(1)*5) + 5;

p = rand(1,l) * 26 - 13;

x = 1:10;
y = ones(size(x));

vals = polyvalnm(p,x,y);
vals_check = polyval(p,x);

pass = check_calc(x,p,vals,vals_check);
if pass ~= 1
    keyboard
end

% case 1b
disp('Case 1b: Polynomial coefficients are single row vector; numel(x) = 1, numel(y) > 1')

l = round(rand(1)*5) + 5;

p = rand(1,l) * 26 - 13;

y = 1:5;
x = ones(size(y));

vals = polyvalnm(p,x,y);
vals_check = polyval(p,x);

pass = check_calc(y,p,vals,vals_check);
if pass ~= 1
    keyboard
end


%% Case 2: coefficients are single column vector 
disp('Case 2a: Polynomial coefficients are single column vector; numel(x) = 1, numel(y) > 1')

l = round(rand(1)*5) + 5;

p = rand(l,1) * 25;

y = 1:10;
x = ones(size(y));

vals = polyvalnm(p,x,y);
vals_check = polyval(p,y);

pass = check_calc(x,p,vals,vals_check);
if pass ~= 1
    keyboard
end

%case 2b
disp('Case 2b: Polynomial coefficients are single column vector; numel(x) > 1, numel(y) = 1')

l = round(rand(1)*5) + 5;

p = rand(l,1) * 25;

x = 1:5;
y = ones(size(x));

vals = polyvalnm(p,x,y);
vals_check = polyval(p,y);

pass = check_calc(y,p,vals,vals_check);
if pass ~= 1
    keyboard
end

%% Case 3. Coeficients x and y non singular dimensions
disp('Case 3a: Polynomial coefficients are non singluar array; numel(x) > 1, numel(y) = 2')

l = round(rand(1)*5) + 5;
% n = round(rand(1)*5) + 5;

p = [zeros(1, l); rand(1,l) * 25];

x = 1:5;
y = ones(size(x));

vals = polyvalnm(p,x,y);

vals_check = polyval(p(2,:),x) .* y;

pass = check_calc(y,p,vals,vals_check);
if pass ~= 1
    keyboard
end

%% Case 4. Coeficients x and y non singular dimensions
disp('Case 4: Polynomial coefficients are non singluar array; numel(x) = 4 (order = 3), numel(y) = 6 (order = 5), fixed size')

n = 5;
m = 3;
p = polyvalnm_coef2mat(rand(polyvalnm_ncoef(n,m),1)*25,n,m);

x = 1:5;
y = ones(size(x));

vals = polyvalnm(p,x,y);

%calculate polynomial suing fixed code. 
%extract polynomial coefficients from p

p00 = p(6,4);
p01  = p(5,4);
p02  = p(4,4);
p03  = p(3,4);
p04  = p(2,4);
p05  = p(1,4);
p10  = p(6,3);
p11  = p(5,3);
p12  = p(4,3);
p13  = p(3,3);
p14  = p(2,3);
p20  = p(6,2);
p21  = p(5,2);
p22  = p(4,2);
p23  = p(3,2);
p30  = p(6,1);
p31  = p(5,1);
p32  = p(4,1);

% calculate polnomial 
vals_check = p00 + p10.*x + p01.*y + p20.*x.^2 + p11.*x.*y + p02.*y.^2 + p30.*x.^3 + p21.*x.^2.*y + p12.*x.*y.^2 + p03.*y.^3 + p31.*x.^3.*y + p22.*x.^2.*y.^2 + p13.*x.*y.^3 + p04.*y.^4 + p32.*x.^3.*y.^2 + p23.*x.^2.*y.^3 + p14.*x.*y.^4 + p05.*y.^5 ;

pass = check_calc(x,p,vals,vals_check);
if pass ~= 1
    keyboard
end

%% POLYVALNM_NCOEF 
% fprintf('\nPOLYVALNM_Ncoef  \n')
% fprintf('Calculate polynomials using Polyvalnm and compare to known functions.\n\n')




fprintf('Solve polynomials using Polyvalnm_solve and compare to known solutions.\n\n')

%% Case 1: coefficients are single row vector 
disp('Polyvalnm_solve, Case 1a: Polynomial coefficients are single row vector; numel(x) > 1, numel(y) = 1')

l = round(rand(1)*5) + 5;

p = rand(1,l) * 26 - 13;

x = 1:10;
% y = ones(size(x));

y_calc = x.^6+x.^5+x.^4+x.^3+x.^2;

vals = polyvalnm_solve(ones(1,6),0,0);
vals_check = polyval(p,x);

pass = check_calc(x,p,vals,vals_check);
if pass ~= 1
    keyboard
end



%% Confirm all tests passed.
fprintf('\nAll test cases passed.\n\n')



end %Polyvalnm_testcases



function out = check_calc(X,coef,val, Val_Check)
out = 1;
fprintf('  Checking sizes of output:')
if isequal(size(Val_Check), size(X))
    fprintf('     Passed \n')
else
    fprintf('     Failed \n')
    out = 0;
end

fprintf('  Checking values of output:')
% if isequal(Val_Check, val)   
if max(sqrt(Val_Check - val).^2 ) < 1E-6 %set like this so the code ignores small rounding errors which may occur in polyvalnm
    fprintf('    Passed \n')
else
    fprintf('    Failed \n')
    out = 0;
end

end %check_calc


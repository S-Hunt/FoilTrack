% Makes a symbolic logic string of the polynomial p with order n, m.
% Where n and m are the orders of the x and y polnomials.
% 
%  syntax:
%       poly_str = polyvalnm_symbolic(p,n,m)


function poly_str = polyvalnm_symbolic(p,n,m)


%check the size of p is correct for the order of the polynomial.
c = polyvalnm_ncoef(n,m);

if c~=sum(~isnan(p(:)))
    error('The order of the polynomial does not correspond to the number of coefficients')
end

p = polyvalnm_coef2mat(p,n,m);
siz = size(p);

p_str = [];
for x = 1:numel(p)
    
    if ~isnan(p(x))
    
        [I,J] = ind2sub(siz,x);
        
        p_term = ['p(',num2str(x),')'];
        if siz(1)-I ~= 0
            p_term_x = ['x^',num2str(siz(1)-I)];
            p_term = [p_term,'*',p_term_x];
        end
        if siz(2)-J ~= 0
            p_term_y = ['y^',num2str(siz(2)-J)];
            p_term = [p_term,'*',p_term_y];
        end
%         x
%         p_term
        
%        p_term = ['p(',num2str(x),')*x^',num2str(siz(1)-I),'*y^',num2str(siz(2)-J)]

        if ~isempty(p_str)
            p_str = [p_str, ' + ', p_term];
        else
            p_str = [p_term];
        end
    end
    
    
    
end

p_str = strrep(p_str, '^1', '');


poly_str = p_str;
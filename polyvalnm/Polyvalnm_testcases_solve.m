%polyvalnm_solve testcases.



%case 1
if 1
    close all
    n = 2;
    m = 1;
    o = polyvalnm_ncoef(n,m);
    p=polyvalnm_coef2mat(ones(1,o),n,m);
    p(2:3) = 0;
    p(5) = 2;
    
    step = 0.1;
    x_range = 0:step:10;
    y_range = -10:step:10;
    [a,b] = meshgrid(y_range, x_range);
    
    s = polyvalnm(p,b,a);
    
    figure
    surf(b,a,s);%, 'EdgeColor', 'none');
    xlabel('y')
    ylabel('x')
    zlabel('I')
    
    
    naughts = polyvalnm_solve(p, 0, x_range);
    
    [~,est_naughts] = min(abs(s),[],2);
    est_naughts = y_range(est_naughts);
    
    figure
    plot(x_range, naughts, '.-', x_range, est_naughts, 'x-k')
    xlabel('y')
    ylabel('location of minimum in x')
    legend('Solution to polynomial', 'estimated from surface')
    
    pause
end

%case2
if 1
    close all
    n = 5;
    m = 2;
    o = polyvalnm_ncoef(n,m);
    p=polyvalnm_coef2mat(ones(1,o),n,m);

    step = 0.05;
    x_range = 0:step:10;
    y_range = -10:step:10;
    [a,b] = meshgrid(y_range, x_range);
    
    s = polyvalnm(p,b,a);
    
    figure
    surf(b,a,s);%, 'EdgeColor', 'none');
    xlabel('y')
    ylabel('x')
    zlabel('I')
    
    
    naughts = polyvalnm_solve(p, 0, x_range);
    
    [~,est_naughts] = min(abs(s),[],2);
    est_naughts = y_range(est_naughts);
    
    
    figure
    plot(x_range, naughts, '.-', x_range, est_naughts, 'x-k')
    xlabel('y')
    ylabel('location of minimum in x')
    legend('Solution to polynomial', 'estimated from surface')
    
    keyboard
end


return
close all
p=polyvalnm_coef2mat(ones(1,9)*1.2345,2,3);
[a,b] = meshgrid(-10:.1:10, -10:.1:10);

figure
s = polyvalnm(p,a,b);
 
 surf(a,b,s)%, 'EdgeColor', 'none');
 
 
 naughts = polyvalnm_solve(p, 0, -10:.1:10);
 
 figure
 plot(-10:.1:10, naughts, '.-')
 
 
function [ F, x, it, time, flag ] = SORl(A, b, x, y, omega, tol)
% This is Exact SOR-like method for Absolute Value Equations 
% (AVE): Ax - |x| = b
% Reference: Ke, Y.F., Ma, C.F.: SOR-like iteration method for solving absolute value equations. Appl. Math. Comput. 311, 195-202 (2017)
%--------------------- Input arguments ------------------------------------

% omega: iteration parameter 

% A, b: data of the problem Ax - |x| = b; 
%       A is nxn real matrix and b is nx1 real vector 
% x,y: initial point; nx1 real vector   
% tol: accuracy for solution of AVE. The method stops at x if || Ax - |x| - b || <= tol
%--------------------- Output arguments -----------------------------------
% x: final iterate
% F: norm of AVE at the final iterate x; 
%    F = || Ax - |x| - b ||
% it: number of iterations of Exact Douglas-Rachford splitting method 
% time: Total CPU time in seconds
% flag = -1: an errror occurred duing the execution
% flag =  0: solution was found
% flag =  1: maximum of iterations reached
% Print initial information
disp('-------------------------------------------------------------------')
disp('------------------   Initial Information   ------------------------')
disp('--------------  SOR-like iteration method -------------------')
disp('-------------------------------------------------------------------')
disp(' ')
disp(' It   || Ax - |x| - b ||')
disp(' ')
tic;

omegam = 1 - omega;
it = 1;
maxit = 1000;%  
FactA = linfactor(A);

while it <= maxit   
    e = y+b;
    res = b + abs(x) - A*x;
    F = norm(res);                   
     if F <= tol
        break
    end
    % Print information
    fprintf(' %3d       %6.2e \n', it, F)
    
    z = linfactor(FactA,e);  
    x1 = omegam*x + omega*z;
    y1 = omegam*y + omega*abs(x1);
    it = it + 1;  
    x = x1;
    y = y1;

end

if it == ( maxit + 1 )
    
    flag = 1;
    time = toc;
    disp('SOR-like iteration method fails to find the solution')
    disp('Maximum of iterations reached ')
    fprintf('Time = %10.4f \n',time)
    
else
    
    flag = 0;
    time = toc;
    disp(' ')
    disp('Solution is found')
    fprintf('Time = %10.4f \n',time)
    
end

end
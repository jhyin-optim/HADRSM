function [F, x, it, time, flag] = Newton(A, b, n, sp, x, tol)

% This is a generalized Newton method for absolute value equations (AVE)
% (AVE): Ax - |x| = b

% Reference: Mangasarian, O.L.: A generalized Newton method for absolute value equations. Optim. Lett. 3(1), 101-108 (2009)

%--------------------- Input arguments ------------------------------------

% n: dimension of the problem

% A, b: data of the problem Ax - |x| = b; 
%       A is nxn real matrix and b is nx1 real vector 

% x: initial point; nx1 real vector 

% sp: sparsity indicator
%     if sp  = 0: full storage organization is used
%     if sp ~= 0: sparse storage organization is used

% tol: accuracy for solution of AVE
%      The method stops at x if ||Ax - |x| - b|| <= tol


%--------------------- Output arguments -----------------------------------

% x: final iterate

% F: norm of AVE at the final iterate x; 
%    F = || Ax - |x| - b ||

% it: number of iterations of Exact Newton method

% time: Total CPU time in seconds

% flag = -1: an error occurred during the execution
% flag =  0: solution was found
% flag =  1: maximum of iterations reached

% Print initial information

disp('-------------------------------------------------------------------')
disp('------------------   Initial Information   ------------------------')
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('--------------  The generalized Newton method  ------------------')
disp('-------------------------------------------------------------------')
disp(' ')
disp(' It   || Ax - |x| - b ||')
disp(' ')

tic;

maxit = 1000;
it = 1;
% 或者改为1
while it <= maxit
    res = b + abs(x) - A*x;
    F = norm(res);
    if F <= tol
        break
    end
    % Print information
    fprintf(' %3d       %6.2e \n',it, F)
    
    if ( sp == 0 )
        AD = A - diag(sign(x));
    else    
        AD = A - sparse(1:n,1:n,sign(x),n,n);
    end
    barF = linfactor(AD); 
    x1 = linfactor(barF,b); 
    it = it + 1;
    x = x1;
end

if it == ( maxit + 1 )
    
    flag = 1;
    time = toc;
    disp('The generalized Newton method fails to find the solution')
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


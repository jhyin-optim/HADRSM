function [F,x,it,time,flag] = DRs(A, b, x, gamma, tol )
% This is the exact Douglas-Rachford splitting method for absolute value equations (AVE) 
% (AVE): Ax - |x| = b
%------------ Input arguments   ---------------------------------
% gamma: iteration parameter   
% A, b: data of the problem Ax - |x| = b; 
%    A is nxn real matrix and b is nx1 real vector 
% x: initial point; nx1 real vector  
% tol: accuracy for solution of AVE  
%   The method stops at x if || Ax - |x| - b || <= tol 
% Output arguments 
% x: final iterate 
% F: norm of AVE at the final iterate x;  
%    F = || Ax - |x| - b || 
% it: number of iterations of Exact Douglas-Rachford splitting method  
% time: Total CPU time in seconds  
% flag = -1: an error occurred during the execution 
% flag =  0: solution was found  
% flag =  1: the maximum of iterations reached 
% Print initial information  
disp('-------------------------------------------------------------------')
disp('------------------   Initial Information   ------------------------')
disp('--------------  Exact Douglas-Rachford splitting method -----------')
disp(' ')
disp(' It   || Ax - |x| - b ||')
disp(' ')
 
tic;  
gamma = 0.5*gamma; 
it = 1; 
maxit = 1000; 
FactA = linfactor(A);

while it <= maxit   
    res = b + abs(x) - A*x;
    F = norm(res);
    if F <= tol
        % Print information
        fprintf(' %3d       %6.2e \n', it, F)
        break
    end
    bar = linfactor(FactA,res);
    x1 = x + gamma*bar;
    it = it + 1;  
    x = x1;
end
if it == ( maxit + 1 )
    
    flag = 1;% flag =  1: the maximum of iterations reached 
    time = toc;
    disp('Exact Douglas-Rachford splitting method fails to find the solution')
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


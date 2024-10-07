function [ F, x, it, time, flag] = IRDRs(A, b, model, x, xs, gamma, tol )
% This is a hybrid acceleration Douglas-Rachford splitting method for absolute value equations (AVE)
% (AVE): Ax - |x| = b
%--------------------- Input arguments   ------------------------------------
% gamma: iteration parameter   
% A, b: data of the problem Ax - |x| = b; 
% A is nxn real matrix and b is nx1 real vector 
% x: initial point; nx1 real vector  
% tol: accuracy for solution of AVE  
%      The method stops at x if || Ax - |x| - b || <= tol 
%--------------------- Output arguments ----------------------------------
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
disp('-------------------------------------------------------------------')
disp('-------------------------------------------------------------------')
disp('---- Hybrid acceleration Douglas-Rachford splitting method -----')
disp('-------------------------------------------------------------------')
disp(' ')
disp(' It   || Ax - |x| - b ||')
disp(' ')
tic;  
alpha = 0.03;  % 0.0132
gamma = 0.5*gamma; 
FactA = linfactor(A);

it = 1; 
maxit = 1000; 
while it <= maxit   
    res = b + abs(x) - A*x;  
    F = norm(res);
    if F <= tol   
        % Print information
        fprintf(' %3d       %6.2e \n', it, F)  
        break
    end
    normsk = norm(x-xs);
    if normsk==0
        alphak = 0;
    elseif model==1
        alphak = min(alpha,1/(it^2*normsk));
    elseif model==2
        alphak = -min(alpha,1/(it^2*normsk));
    elseif model==3
        alphak = (-1)^it*min(alpha,1/(it^2*normsk));
    else
        alphak = alpha;
    end
    y = x + alphak*(x-xs); 
    res1 = b + abs(y) - A*y;
    F = norm(res1);
    if F <= tol   
        % Print information
        fprintf(' %3d       %6.2e \n', it, F) 
        break
    end
    bary = linfactor(FactA,res1); 
    x1 = y+gamma*bary;  
    it = it + 1;
    xs = x;
    x = x1;
end
if it == ( maxit + 1 )
    flag = 1; % flag =  1: the maximum of iterations reached 
    time = toc;
    disp('Hybrid acceleration Douglas-Rachford splitting method fails to find the solution')
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
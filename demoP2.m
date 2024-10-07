clc;
close all

% set random number seed
rng(2016)

% setup TXT document
fid_tex=fopen('mytextP2.txt','w'); 
problem_set = [100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000]; % 5000 6000 7000 8000];
np = length(problem_set);

% parameter
gamma = 1.91;
omega = 0.94;
tol   = 1e-7;

% run
for index=1:np
    n = problem_set(index);
    x_star = (-1).^(1:n)'; 
    progress_r = [];
    for repeats = 1:5
        V = 2*rand(n,1)+1;
        V = diag(V);
        W = randn(n);
        U = orth(W);
        A = (U*V)*U';
        bvector = A*x_star - abs(x_star);
        % input the initial point
        x0 = randn(n,1); %rand(n,1)-rand(n,1); %-100+200*rand(n,1);
%         normAinv = norm(inv(A));

        [Fd, xd, itd, timed, ~] = DRs(A, bvector, x0, gamma, tol);  
    
        [Fs, xs, its, times, ~] = SORl(A, bvector, x0, x0, omega, tol);
    
%         [Fn, xn, itn, timen, ~] = Newton(A, bvector, n, 1, x0, tol);
    
        [Fp, xp, itp, timep, ~] = Picard(A, bvector, tol);
    
        [Fir1, xir1, itir1, timeir1, ~] = IRDRs(A, bvector, 1, x0, x0, gamma, tol);
    
        [Fir2, xir2, itir2, timeir2, ~] = IRDRs(A, bvector, 2, x0, x0, gamma, tol);
    
        [Fir3, xir3, itir3, timeir3, ~] = IRDRs(A, bvector, 3, x0, x0, gamma, tol);
        
        progress_r = [progress_r;itir1,timeir1,Fir1,itir2,timeir2,Fir2,itir3,timeir3,Fir3,itd,timed,Fd,its,times,Fs,itp,timep,Fp];%itn,timen,Fn,];
    end
    TM = mean(progress_r); 
    fprintf(fid_tex,'%d & %.1f/%.4f/%.2e & %.1f/%.4f/%.2e & %.1f/%.4f/%.2e & %.1f/%.4f/%.2e\n & %.1f/%.4f/%.2e & %.1f/%.4f/%.2e\\\\ \r\n', ... 
                n,TM);
end
%% ¹Ø±ÕÎÄ¼þ
fclose(fid_tex);
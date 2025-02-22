clc;
close all

% set random number seed
rng(2016)

% setup TXT document
fid_tex=fopen('mytextP4_04.txt','w'); 
problem_set = [100 200 300 400 500 600 700 800 900 1000 2000 3000];
np = length(problem_set);

% parameter
gamma = 1.93;
omega = 0.94;
tol   = 1e-7;

% run
for index=1:np
    n = problem_set(index);
    progress_r = [];
    for repeats = 1:50
        x_star = 200*(rand(n,1)-0.5);
        rc = rand(n,1);
        rc = rc./(min(rc)*rand);
        A = sprand(n,n,0.4,rc);
        bvector = A*x_star - abs(x_star);
        % input the initial point
        x0 = 200*(rand(n,1)-0.5); 

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
%% 关闭文件
fclose(fid_tex);
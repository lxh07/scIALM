mask = csvread("preprocess/data/PBMC/masked/PBMC_01.csv",1,1);

newDataFlag = 1;
if newDataFlag == 1
    
    clc ;
    close all ;
    %m = 2843 ;
    pdr = 6;
    %n = 13003 ;
    [m, n]=size(mask);
    r = 10;
    [i,j] = find(mask);
    p = length(i);   
    rho_s = p / (m * n);
    
    [I J col omega] = allsample(m, n, p, mask);
    V = UVtOmega(mask, I, J, col);
    
    D = spconvert([I,J,V; m,n,0]);
%     clear I J col;
end

disp('m');
disp(m);
disp('n');
disp(n);
disp('r');
disp(r);
disp('p');
disp(p);
disp('pdr');
disp(pdr);
disp('rho_s');
disp(rho_s);

tic;
[A iter svp] = inexact_alm_mc(D, 1e-4, 1000);
tElapsed = toc;

imputed = A.U*A.V';
data_table = array2table(imputed);disp('Time');
disp(tElapsed);
writetable(data_table, "preprocess/data/PBMC/result/PBMC_01.csv");

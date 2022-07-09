% Global sensitivity analysis
% Model output is the kinectics of the tumor proliferating cells + damaging cells(100 mg/kg weekly )
% step1. choose the most influential parameter
% step2. Variance-based sensitivity analysis
% purturb the parameter with the range (determine the range)
% create 1000 parameter sets using Latin hypercube sample
close all
clear 
dbstop if error
%use Morries method to select subsets of parameters that seemed likely to have substantial effects on model outputs.
num_sample = 1000; %Total number of samples is 2000
theta_PK =   [329.539000000000,320.071000000000,16823.4720000000,2999.82100000000,17334.9770000000];
theta_cellular =[0.980860000000000;538.450690000000;171.153100000000;536.641620000000;6.70454000000000;135.995300000000;27.1384400000000;2325.15421000000;0.0890600000000000;0.521000000000000;0.114470000000000;0.310230000000000;0.0199000000000000;0.839790000000000;5.53387000000000;0.959080000000000;0.807740000000000;0.554160000000000;0.395870000000000;4.35978000000000;0.0385200000000000;2.88205000000000;0.152570000000000;21.9351600000000;0.154750000000000;17.4304600000000;0.0321300000000000];
theta_DSB = [0.212940000000000;32.4999300000000;3.71261000000000;73.3892700000000;3.69472000000000;0.303530000000000;3.33921000000000;1.38163000000000;0.276720000000000;35.6500000000000;6.03853000000000;0.0793900000000000;1.80119000000000;10.5164600000000;151.918160000000;10.0816500000000;194.020220000000;0.121190000000000;2.80318000000000];
theta_TGI = [497.368268819054;7.63003855255054;0.673708616383969;5.13957298019130;7.67657941208417e-05;4.84252683969523;1.16636621985418e-05];
num_para = length(theta_PK) +   length(theta_cellular)  + length(theta_DSB)  + length(theta_TGI);
A = lhsdesign(num_sample,num_para);
B = lhsdesign(num_sample,num_para);
day2min = 24*60;
T = (5:5:45)*day2min;
T_num = length(T);
mu = [theta_PK(:)   ;  theta_cellular(:) ;    theta_DSB(:) ; theta_TGI(:) ;  ];
mu = mu';
ub_col = zeros(1, num_para);
lb_col = zeros(1, num_para);
for m = 1:num_para
      ub_col(m) = mu(m)*1.5; 
      lb_col(m) = mu(m)*0.5;
end
A = lhsdesign(num_sample,num_para);
B = lhsdesign(num_sample,num_para);
L_A =  repmat(lb_col, num_sample,1)  +  A.*repmat(ub_col- lb_col, num_sample,1 )  ;
L_B =  repmat(lb_col, num_sample,1)  +  B.*repmat(ub_col- lb_col, num_sample,1 )  ;
%loop over the parameters
Si_col = zeros(num_para,T_num);
[Tsim, Ysim,f_A_prime] =  kinetics_TGI(mu,1,T);
%h = parallelcoords(L_A(:,1:5) )
for i = 1: num_para  %loop over the parameters
    M_A = zeros(num_sample*2,T_num);
    M_B = zeros(num_sample*2,T_num);
    M_A_prime = zeros(num_sample*2,T_num);
    for j = 1:num_sample
        theta_A = L_A(j,:);
        theta_B = L_B(j,:);
        [~,~,f_A] =  kinetics_TGI(theta_A,1,T);
        [~,~,f_B] =  kinetics_TGI(theta_B,1,T);
        if i == 1
            theta_A_prime = [ theta_B(1)   theta_A(2:end)];
        else
            theta_A_prime = [theta_A(1:i-1)  theta_B(i)   theta_A(i+1:end)];
        end
        [Tsim, Ysim,f_A_prime] =  kinetics_TGI(theta_A_prime,1,T);
        M_A(j,:) =   f_A;
        M_B(j,:) =   f_B;
        M_A_prime(j,:) = f_A_prime;
        if rem(j,20) == 0 %100*floor(curr_T/100) = curr_T
            disp(j)
        end
    end
end
for t = 1:length(T)
    temp =  M_A(:,t).*(M_A_prime(:,t)- M_B(:,t)) ;
    Si = mean(temp);
end
    

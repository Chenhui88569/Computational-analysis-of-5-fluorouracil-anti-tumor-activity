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
% theta_PK =   [329.539000000000,320.071000000000,16823.4720000000,2999.82100000000,17334.9770000000];
% theta_cellular =[0.980860000000000;538.450690000000;171.153100000000;536.641620000000;6.70454000000000;135.995300000000;27.1384400000000;2325.15421000000;0.0890600000000000;0.521000000000000;0.114470000000000;0.310230000000000;0.0199000000000000;0.839790000000000;5.53387000000000;0.959080000000000;0.807740000000000;0.554160000000000;0.395870000000000;4.35978000000000;0.0385200000000000;2.88205000000000;0.152570000000000;21.9351600000000;0.154750000000000;17.4304600000000;0.0321300000000000];
% theta_DSB = [0.212940000000000;32.4999300000000;3.71261000000000;73.3892700000000;3.69472000000000;0.303530000000000;3.33921000000000;1.38163000000000;0.276720000000000;35.6500000000000;6.03853000000000;0.0793900000000000;1.80119000000000;10.5164600000000;151.918160000000;10.0816500000000;194.020220000000;0.121190000000000;2.80318000000000];
% theta_TGI = [497.368268819054;7.63003855255054;0.673708616383969;5.13957298019130;7.67657941208417e-05;4.84252683969523;1.16636621985418e-05];
% num_para = length(theta_PK) +   length(theta_cellular)  + length(theta_DSB)  + length(theta_TGI);
num_worker = 190;
num_sample = 500; %Total number of samples is 500
load('Para.mat','Para_Info');
Para_Name = Para_Info.Para_Name;
Para_Value =  Para_Info.Para_Value;
num_para  = length(Para_Name);
para_change = ["K_m";"K_{ m , efflux }";"K_{ m, 54}";"K_{m, influx }";"P_{max}";"V_{ max ,  efflux }";"V_{ max , 54}";"V_{max , 75}";"V_{max , influx }";"V_{max}";"\lambda_{d}";"\lambda_{g}";"k_{03}";"k_{56}";"k_{59}";"k_{65}";"k_{95}"];
num_para_change = length(para_change);
para_change_idx = zeros(num_para_change,1);
for i = 1:num_para_change
    para_change_idx(i) = find(strcmp(Para_Name ,para_change(i))); %find indices of the parameters to be sampled in the index para_change.
end
A = lhsdesign(num_sample,num_para_change);
B = lhsdesign(num_sample,num_para_change);
day2min = 24*60;
T = (5:5:45)*day2min;
T_num = length(T);
ub_col = zeros(1, num_para_change);
lb_col = zeros(1, num_para_change);
mu_para_change = Para_Value(para_change_idx);
for m = 1: num_para_change
      ub_col(m) = mu_para_change(m)*2; 
      lb_col(m) = mu_para_change(m)*0.5;
end
L_A =  repmat(lb_col, num_sample,1)  +  A.*repmat(ub_col- lb_col, num_sample,1 )  ;
L_B =  repmat(lb_col, num_sample,1)  +  B.*repmat(ub_col- lb_col, num_sample,1 )  ;
%% calculate the overall response variance 

%loop over the parameters
Si_col = zeros(num_para_change,T_num);
%[Tsim, Ysim,f_A_prime] =  kinetics_TGI(mu,1,T);

%h = parallelcoords(L_A(:,1:5) )
M_A = zeros(num_sample,T_num); % matrix A 
M_B = zeros(num_sample,T_num); % matrix B

parpool(num_worker)
parfor j = 1:num_sample
    theta_A = L_A(j,:);
    theta_A_complete = Para_Value;
    theta_A_complete_var = theta_A_complete ;   theta_A_complete_var(para_change_idx) = theta_A;
    theta_B = L_B(j,:);
    theta_B_complete = Para_Value;
    theta_B_complete_var = theta_B_complete ;   theta_B_complete_var(para_change_idx) = theta_B;
    [~,~,f_A] =  kinetics_TGI(theta_A_complete_var ,1,T);
    [~,~,f_B] =  kinetics_TGI(theta_B_complete_var ,1,T);
    M_A(j,:) =   f_A;
    M_B(j,:) =   f_B;
end
delete(gcp('nocreate'))

%% overall response variance
M_col  =[M_A;
    M_B];
V_col = mean(M_col.^2,1) - mean(M_col,1).^2;
V_col = abs(V_col );
%%

M_A_prime_col = cell( num_sample ,1);
parpool(num_worker)
parfor j = 1 : num_sample  %loop over the parameters
    M_A_prime = zeros(num_para_change, T_num);
    for i = 1: num_para_change
        % To prevent the error 'unclassified variable.', create a temporary array in each paraloop and then fill in in the nested for loop.

        theta_A = L_A(j,:);
        theta_A_complete = Para_Value;
        theta_A_complete_var = theta_A_complete ;   theta_A_complete_var(para_change_idx) = theta_A;
        theta_B = L_B(j,:);
        theta_B_complete = Para_Value;
        theta_B_complete_var = theta_B_complete ;   theta_B_complete_var(para_change_idx) = theta_B;
        para_id = para_change_idx(i);
        if  para_id == 1
            theta_A_prime = [  theta_B_complete_var(1)  ;  theta_A_complete_var(2:end)];
        else
            theta_A_prime = [theta_A_complete_var(1: para_id-1); theta_B_complete_var( para_id) ;  theta_A_complete_var( para_id+1:end)];
        end
        [Tsim, Ysim,f_A_prime] =  kinetics_TGI(theta_A_prime,1,T);
        M_A_prime(i,:) = f_A_prime;
    end
    M_A_prime_col{j} = M_A_prime;
end
delete(gcp('nocreate'))

%% calculate the overall response variance 
Si_col = zeros(num_para_change,T_num);
STi_col = zeros(num_para_change,T_num);
for indi_para = 1:num_para_change
    M_A_prime = zeros(num_sample,T_num);
    for m = 1:num_sample
        temp =  M_A_prime_col{m};
        M_A_prime(m,:) = temp( indi_para,:);
    end
    for t = 1:length(T)
        temp_i =  M_B(:,t).*(M_A_prime(:,t)- M_A(:,t)) ;
        Si = mean(temp_i);
        temp_Ti = (M_A_prime(:,t)- M_A(:,t)).^2;
        STi  = 1/2*mean(temp_Ti);
        Si_col(indi_para,t) = Si;
        STi_col(indi_para,t) =  STi; 

    end
end

newdir = 'Sobol';
if ~isfolder(newdir)
    mkdir(newdir);
else
    rmdir(newdir,'s');
    mkdir(newdir)
end

f_s = figure;
set(f_s, 'Position', get(0, 'Screensize'));
T_trun_idx = [1, ceil(9/2),9];
STi_col= abs(STi_col); Si_col = abs(Si_col);
for indi_T = 1:length(T)
    Si_col_indiT =  Si_col(:,indi_T)/V_col(indi_T) ;
    STi_col_indiT = STi_col(:,indi_T)/V_col(indi_T)  ;
    if ismember(indi_T, T_trun_idx )
        subplot(3,1,find(T_trun_idx == indi_T))
        bar( [Si_col_indiT  STi_col_indiT ])
        ylabel('Sobol Index')
        xlabel('Sensitivity input')
        xticks(1:num_para_change)
        xticklabels(para_change)
        legend({'First order','Total order'},'Location','northwestoutside')
        title(['Sobol index with tumor volume at ' num2str(T(indi_T)/day2min) ' th day as output'     ])
    end
end
saveas(f_s, fullfile(newdir , 'Sobol.tif'))
saveas(f_s, fullfile(newdir , 'Sobol.fig'))
save(fullfile(newdir, 'Sobol_info.mat'), 'M_col',  'M_A_prime_col') 

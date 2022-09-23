clear
close all
t_plasma  = [0;10;20;30;60;90;120;240];
T= 10000;
%  V_max = theta_plasma(1);
%     Q_21 = theta_plasma(2);
%     V_1	= theta_plasma(3);
%     V_2	= theta_plasma(4);
%     K_m = theta_plasma(5);
%StepSize = [1000 5 100 500  1000  10e4 ];


init_para = [   9.640e-1
5.400e+2
1.710e+2
5.370e+2
6.350e+0
1.270e+2
2.520e+1
2.320e+3
 1.300e-2
4.820e-1
6.816e-2
6.130e-4
1.737e-2
8.690e-1
5.786e+0
9.615e-1
9.614e-1
5.881e-1
7.391e-1
4.540e+0
3.913e-2
2.886e+0
1.753e-1
2.194e+1
1.984e-1
1.740e+1
3.234e-2
2e2
2e2
5
2
0.2
2
];
init_para = init_para';
num_para = length(init_para);

para_label = [ "Q_{31}"
"V_{max , influx }"
"K_{m, influx }"
"V_{ max ,  efflux }"
"K_{ m , efflux }"
"k_{03}"
"V_{ max , 54}"
"K_{ m, 54}" 
"k_{65}"
"k_{56}"
"k_{06}"
"\gamma_{lag, RNA }"
"V_{max , 75}"
"K_{m, 75}"
"V_{max , 57}"
"K_{m, 57}"
"k_{07}"
"\gamma_{lag, DNA}"
"T_{d,RNA}"
"T_{d,DNA}"
"k_{95}"
"k_{59}"
"k_{09}"
"K_{dUMP}"
"G_{0}"
"k_{cat}"
"k_{08}"
"\sigma_{intra}"
"\sigma_{ana}"
"\sigma_{RNA}"
"\sigma_{DNA}"
"\sigma_{TS}"
"\sigma_{dUMP}"
];
% 
load('Para_Cellular_col_uni_new.mat')
Para_col_hist = Para_col;
Para_col = zeros( T,    num_para ) ; 
Para_col(1,:) =  Para_col_hist(end,1:num_para); % init_para;%Para_col_hist(end,1:num_para); % init_para ; % Para_col_hist(end,1:num_para);
% 
% S_curr = zeros(num_para);
% for i = 1: num_para
%     S_curr(i,1:i) =  abs(randn);
%     if i > 8
%         S_curr(i,1:i) =  abs(randn)*0.01;
%     end
% end

S_curr  = S_curr( 1: num_para, 1:num_para );
myfun = @Target_fun_LE_Cellular_Uniform;
[Para_col,  acc_rate_col , S_curr] = MyRAM(  T, num_para, Para_col , S_curr,myfun,para_label,[]);
Para_col = [Para_col_hist; Para_col(2:end,:)];
outputfile_dir = 'MCMC_cellular/';
save( [ outputfile_dir   'Para_Cellular_col_uni_new_2.mat']  ,'Para_col'  , 'acc_rate_col',  'S_curr' ,'-v7.3');
    
   
clear
close all
T =  5000;

init_para = [   
 2.135e-1
 3.250e+1
 3.714e+0
 7.334e+1
 3.770e+0
 3.026e-1
 3.335e+0
 1.353e+0
 2.761e-1
 3.567e+1
 6.045e+0
 7.853e-2
 1.790e+0
 1.051e+1
 1.519e+2
 1.007e+1
 1.940e+2
 1.255e-1
 2.794e+0
 0.05
 2e2
];
init_para = init_para';
num_para = length(init_para);

para_label = [ "k_1"
"k_2"
"k_3"
"k_4"
"k_5"
"k_6"
"k_7"
"k_8" 
"k_9"
"k_B"
"k_A"
"k_{10}"
"T_{d,DSB}"
"V_{max,dNTP}"
"K_{m,dNTP}"
"V_{max,HR}"
"K_{m,HR}"
"k_i"
"k_0"
"\sigma_{dNTP}"
"\sigma_{DSB}"

];

load('Para_dNTP_DSB_col_2.mat')
Para_col_hist = Para_col;
Para_col = zeros( T,    num_para ) ; 
Para_col(1,:) =  Para_col_hist(6496,:); %init_para ;% Para_col_hist(end,1:num_para);

% S_curr = zeros(num_para);
% for i = 1: num_para
%     S_curr(i,1:i) =  abs(randn)*0.01;
% end

S_curr  = S_curr( 1: num_para, 1:num_para );
myfun = @Target_fun_LE_dNTP_DSB;
[Para_col,  acc_rate_col , S_curr] = MyRAM(  T, num_para, Para_col , S_curr,myfun,para_label);
Para_col = [Para_col_hist(1:6495,:); Para_col];
save Para_dNTP_DSB_col_2.mat  Para_col    acc_rate_col  S_curr
    
   
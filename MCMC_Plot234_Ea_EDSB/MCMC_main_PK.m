clc
close all
t_plasma  = [0;10;20;30;60;90;120;240];
T = 10000;
num_para = 5+1;

load Para_PK_col.mat Para_col S_curr
Para_col_hist = Para_col;
Para_col = zeros( T,    num_para ) ; 

Para_col(1,:) = Para_col_hist(end,:);

[Para_col,  acc_rate_col , S_curr] = MyRAM(  T, num_para, Para_col , S_curr);
Para_col = [Para_col_hist; Para_col];
save Para_PK_col.mat  Para_col   S_curr
    
   
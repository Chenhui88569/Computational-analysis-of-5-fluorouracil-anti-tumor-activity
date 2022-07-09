clear
close all
T =6000;

%init_para = [2.300e-4 1.320e+1 5.577e-6,1,1];
init_para = [1300  69.380333945 0.526513177 49.468042148 0.000249874 7.717727070 0.000004283 0.5 0.5 ];
%init_para = init_para';
Model_Index = 2; 
if ismember( Model_Index, [4,6,7])
    para_label = [  "\lambda_{g}",   "P_{max}",   "\lambda_{d}",  "\sigma_{treated}","\sigma_{control}" ];
else
    para_label = [ "T_{d,Tv}", "IC_{50}",   "E_{max , damage }",   "EC_{50,damage}",   "\lambda_{g}",   "P_{max}",   "\lambda_{d}",  "\sigma_{treated}","\sigma_{control}" ];
end
% 
num_para = length(para_label);
load("Para_TGI_col_Model2.mat")
Para_col_hist = Para_col;
Para_col = zeros( T,    num_para ) ; 
Para_col(1,:) = Para_col_hist(1,:);%Para_col_hist(end,1:num_para); %init_para ;% Para_col_hist(end,1:num_para);
Para_col(1,end-1) = 0.5;
Para_col(1,end) = 0.5;
S_curr = zeros(num_para);
for i = 1: num_para
    S_curr(i,1:i) =  abs(randn);
    if  i == 1
          S_curr(i,1:i) =  abs(randn)*4;
    end
    if  i == 5 ||  i == 7
         S_curr(i,1:i) =  abs(randn)*1e-4;
    end
    if  i == 2
        S_curr(i,1:i) =  abs(randn);
    end
end

% S_curr = zeros(num_para);
% for i = 1: num_para
%     S_curr(i,1:i) =  abs(randn)*5;
%     if  i == 1 ||  i == 3
%         S_curr(i,1:i) =  abs(randn)*1e-3;
%     end
% end

S_curr  = S_curr( 1: num_para, 1:num_para );
myfun = @Target_fun_LE_TGI;
[Para_col,  acc_rate_col , S_curr] = MyRAM(  T, num_para, Para_col , S_curr,myfun,para_label,Model_Index );
%Para_col = [Para_col_hist; Para_col];
save Para_TGI_col_Model2_test2.mat  Para_col  acc_rate_col S_curr
  

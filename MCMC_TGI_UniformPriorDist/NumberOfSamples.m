close all
clear
portion_coll = [0.7; 0.3;0.3;0.07;0.2;0;0.15;0;0.15 ];
num_sample_col = ones(1,9);
for Model_Index = 1:9
    load(['MCMC_TGI_' num2str(Model_Index) '/' 'Para_TGI_col_uni_' num2str(Model_Index) '.mat'])

    %load("Para_TG_col_20TwoTimesAWeek_3Weeks.mat")
    Para_col_for_analysis = Para_col;
    
    if ismember( Model_Index, [4,6,7])
        para_label = [  "\lambda_{g}",   "P_{max}",   "\lambda_{d}",  "\sigma_{treated}^2","\sigma_{control}^2" ];
    else
        para_label = [ "T_{d,Tv}", "IC_{50}",   "E_{max , damage }",   "EC_{50,damage}",   "\lambda_{g}",   "P_{max}",   "\lambda_{d}",  "\sigma_{treated}^2","\sigma_{control}^2" ];
    end
    percent = portion_coll(Model_Index);%0.3 96
    num_sample = size(Para_col_for_analysis,1);
    num_para = length(para_label);
    Para_trun = Para_col_for_analysis( floor(num_sample*percent)+1:end,1:num_para);
    Para_mean_col = zeros(num_para,1);
    num_sample_truncated =  size( Para_trun,1);
    num_sample_col(Model_Index) =  num_sample_truncated ;
end
num_sample_col_pluecolon = append(string(num_sample_col),',');
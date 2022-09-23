%This script is used to organize the parameters for nine TGI models in a table.
% T_col containing confidence interval is used in supplementay material 3
% MeanPara_array with only mean parameter values is used in manuscript.
close all
clear
portion_coll = [0.7; 0.3;0.3;0.07;0.2;0;0.15;0;0.15 ];
T_col =  cell(9, 9);
Sample_num_col = zeros(1,9);
MeanPara_array = zeros(9,7);
sigma_col   = zeros(9,2);
MeanPara_CI_array =  cell(7,9);
for Model_Index = 1:9
    load(['MCMC_TGI_' num2str(Model_Index) '/' 'Para_TGI_col_uni_' num2str(Model_Index) '.mat'])

    
    Para_col_for_analysis = Para_col;
    if ismember( Model_Index, [4,6,7])
        para_label = [  "\lambda_{g}",   "P_{max}",   "\lambda_{d}", "\sigma_{treated}","\sigma_{control}" ];
    else
        para_label = [ "T_{d,Tv}", "IC_{50}",   "E_{max , damage }",   "EC_{50,damage}",   "\lambda_{g}",   "P_{max}",   "\lambda_{d}",  "\sigma_{treated}","\sigma_{control}" ];
    end
    percent = portion_coll(Model_Index);%0.3 96
    num_sample = size(Para_col_for_analysis,1);
    num_para = length(para_label);
    Para_trun = Para_col_for_analysis( floor(num_sample*percent)+1:end,1:num_para);
    Para_mean_col = zeros(num_para,1);
    Sample_num_col(Model_Index ) = size( Para_trun ,1);
    for idx = 1: num_para 
        %subplot(2,3, i)
        dbstop if error
        Para_mean = mean(Para_trun(:,idx ));
        %Para_mean = round(Para_mean);
        Para_mean_col(idx ) =  Para_mean;
        x = Para_trun(:,idx );
        CI  = quantile(x, [0.25, 0.75]);
        MeanPara_CI_array{idx,Model_Index}= [  num2str(Para_mean)  '('    num2str(CI(1) ) ','      num2str(CI(2))  ')'     ];
        %bar_label = strcat( '(',para_label(idx ),')','_{avg}' );
        %temp = strcat(num2str( Para_mean) ,'(', num2str(CI(1)), ',', num2str(CI(2)),')');
        % four digits after the decimal point.
       temp1 = sprintf('%0.4fe%i', 10^mod(log10(Para_mean),1),floor(log10(Para_mean)) );
       temp2 = sprintf('(%0.4fe%i,%0.4fe%i)', 10^mod(log10(CI(1)),1),floor(log10(CI(1))),...
            10^mod(log10(CI(2)),1),floor(log10(CI(2))) );
        if idx <=  7 % for the other models
            if  ~ismember( Model_Index, [4,6,7])
                MeanPara_array(Model_Index,idx) = Para_mean;
            end
        end
        if  idx <= 4 % model 4,6,7
                if  ismember( Model_Index, [4,6,7])
                    MeanPara_array(Model_Index, idx ) =  MeanPara_array(3, idx ) ;
                    if idx  <= 3
                        MeanPara_array(Model_Index, idx +4 ) =  Para_mean;
                    end
                end
        end
        %0.2f x 10^%i', 10^mod(log10(surfaceIons),1),floor(log10(surfaceIons))
      if ismember( Model_Index, [4,6,7])
            if  ismember( idx, 1:4)
               T_col{ Model_Index, idx} = 'b';
            end
            T_col{ Model_Index, idx+4} =    [temp1 '\\' temp2];
      else
        T_col{ Model_Index, idx} =[temp1 '\\' temp2];
      end
    end
    sigma_col(Model_Index,:) = Para_mean_col(end-1:end);
end
T_col = table(T_col)
%writetable(T_col, 'AllTGI_para_col.xlsx')
%save NineTGI_Para.mat  MeanPara_array MeanPara_CI_array
  
% For individual TGI
close all
clear
Model_Index = 7;
%load(['Para_TGI_col_Model' num2str(Model_Index) '.mat'])
load(['MCMC_TGI_' num2str(Model_Index) '/' 'Para_TGI_col_uni_' num2str(Model_Index) '.mat'])
portion_coll = [0.7; 0.3;0.3;0.07;0.2;0;0.15;0;0.15 ];
Para_col_for_analysis = Para_col;
% It's used in the situation where the MCMC program is terminated prematurely.
if find(Para_col_for_analysis(:,1) == 0)
    idx_nonzero = find(Para_col_for_analysis(:,1) == 0);
    Para_col_for_analysis =  Para_col_for_analysis( 1:   idx_nonzero-1 ,: );
    Para_col = Para_col_for_analysis;
    save(['MCMC_TGI_' num2str(Model_Index), '/', 'Para_TGI_col_uni_' num2str(Model_Index) '.mat'] , 'Para_col', 'acc_rate_col',  'S_curr' ,'-v7.3' )
end
    
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

factor = 9;
T_col =  cell(num_para,1);
for iFig = 1: ceil(num_para/factor )
    f= figure;
    set(0, 'CurrentFigure', f)
    set(f, 'Position', get(0, 'Screensize'));
    for i = 1:factor 
        subplot(3,3, i)
        idx = factor*(iFig-1)+i;
        dbstop if error
        histogram(Para_trun(:,idx  ) )
        hold on
        [counts, edges] = histcounts(Para_trun(:,idx ));
        locs = movmean(edges, 2, 'Endpoints', 'discard');
        plot(locs, counts, 'LineWidth', 3);
        ax = gca;
        title(para_label(idx ),'FontSize', 15)
        Para_mean = mean(Para_trun(:,idx ));
        %Para_mean = round(Para_mean);
        Para_mean_col(idx ) =  Para_mean;
        hold on
        xline(Para_mean,'LineWidth',3)
        x = Para_trun(:,idx );
        CI  = quantile(x, [0.25, 0.75]);        %CI = round(CI,3);
        %bar_label = strcat( ' bar', '{', para_label(i),'}'   );
        bar_label = strcat( '(',para_label(idx ),')','_{avg}' );
        temp = strcat(num2str( Para_mean) ,'(', num2str(CI(1)), ',', num2str(CI(2)),')');
        T_col{idx} = temp;

        %text( ax.XTick(1),  (ax.YTick(end)+ax.YTick(end-1))/2  ,temp,'tex'  )
        hold off
        if  idx == num_para
            break
        end
    end
end
T = table(T_col);

for  i = 1: floor(num_para/6)+1
       f =  figure;
       set(0, 'CurrentFigure', f);
      % set(f, 'Position', get(0, 'Screensize'));
        for j = 1:6
            subplot(2,3, j)
            idx = 6*(i-1)+j;
            plot(Para_trun(:, idx ),'k-','LineWidth',1.5)
            xlabel('iteration')
            %axis([1+T*fix(curr_T/T),T+T*fix(curr_T/T),0,15]);
            grid off
            title(para_label(idx ),'FontSize', 15)
            if  idx == num_para
                break
            end
        end
end

C_diff_cell =  kinetics_TumorVolume_Treated(Para_mean_col,Model_Index);
 lss = vecnorm([C_diff_cell{1}; C_diff_cell{2}],2)
 Plot_TGIModel(Para_trun,Para_mean_col,Model_Index);
 

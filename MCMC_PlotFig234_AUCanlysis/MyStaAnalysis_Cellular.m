close all
clear
load("Para_Cellular_col.mat")
Para_col_for_analysis = Para_col; %num_sample *num_variable

newdir_parent= 'Plot_Cellular_UpTodUMP_DSB';   % Your destination folder 
newdir  = fullfile(newdir_parent , 'Cellular_UpTodUMP'  );
if ~isfolder( newdir   ) 
   mkdir( newdir   ) ;
else 
  rmdir( newdir ,'s'  ) ;
  mkdir( newdir ) ;
end

para_label = [ "Q_{31}"
"V_{max , influx }"
"K_{m, influx }"
"V_{max, efflux }"
"K_{m, efflux }"
"k_{03}"
"V_{max, 54}"
"K_{m, 54}" 
"k_{65}"
"k_{56}"
"k_{06}"
"\gamma_{lag, RNA }"
"V_{max, 75}"
"K_{m, 75}"
"V_{max , 57}"
"K_{m, 57}"
"k_{07}"
"\gamma_{lag, DNA}"
"T_{d, RNA}"
"T_{d, DNA}"
"k_{95}"
"k_{59}"
"k_{09}"
"K_{dUMP}"
"G_{0}"
"k_{cat}"
"k_{08}"
"\sigma_{intra}"
"\sigma_{FNUC}"
"\sigma_{RNA}"
"\sigma_{DNA}"
"\sigma_{TS}"
"\sigma_{dUMP}"
];
percent = 0.25;%0.26 
num_sample = size(Para_col_for_analysis,1);
num_para = length(para_label);
Para_mean_col = zeros(num_para,1);
Para_trun = Para_col_for_analysis( floor(num_sample*percent)+1:end,1:num_para);
num_trun = size(Para_trun,1);
T_col =  cell(num_para ,2); T_onlyest =  cell(num_para ,2);
T_col(:,2) = num2cell(para_label );
T_onlyest(:,2) = num2cell(para_label );
T_med_col = zeros(num_para,1);
factor =  35;
for iFig = 1: ceil(num_para/factor )
    f= figure;
    set(0, 'CurrentFigure', f)
    set(f, 'Position', get(0, 'Screensize'));
    for i = 1:factor 
        subplot(7,5, i)
        idx = factor *(iFig-1)+i;
        dbstop if error
        %% plotting histogram
        histogram(Para_trun(:,idx  ) )
        hold on
        [counts, edges] = histcounts(Para_trun(:,idx ));
        locs = movmean(edges, 2, 'Endpoints', 'discard');
        plot(locs, counts, 'LineWidth', 3);
        ax = gca;
        title(para_label(idx ),'FontSize', 15)
        Para_mean = mean(Para_trun(:,idx ));
        Para_mean = round(Para_mean,5);
        Para_mean_col(idx ) =  Para_mean;
        hold on
        xline(Para_mean,'LineWidth',3)
        T_onlyest{idx,1} =  num2str(Para_mean);
        x = Para_trun(:,idx );
        %% Calculate credible interval
        CI  = quantile(x, [0.25, 0.75]);
        %bar_label = strcat( ' bar', '{', para_label(i),'}'   );
        %bar_label = strcat( '(',para_label(idx ),')','_{avg}' );
        temp = strcat( num2str( Para_mean) ,'(', num2str(CI(1)), ',', num2str(CI(2)),')');
        T_col{idx,1} = temp;
        T_med_col(idx) = median(x);
        %text( ax.XTick(1),  (ax.YTick(end)+ax.YTick(end-1))/2  ,temp,'tex'  )
        hold off
        if  idx == num_para
            break
        end
    end
    sgtitle( 'Histograms for parameters in \Theta_2','FontSize', 15)
    saveas( f, fullfile(newdir  , ['hist'  num2str(iFig) '.tif']));

end
%T = table(T_col);
%T_onlyest = table(T_onlyest);
%% covariance matrix 
num_sample_afterburnin = size(Para_trun,1);
Cov = zeros(num_para, num_para );
for i = 1:num_sample_afterburnin
    Cov  = Cov + (Para_trun(i,:) - Para_mean_col')' * (Para_trun(i,:) - Para_mean_col');
end
Cov = 1/(num_sample_afterburnin-1)*Cov; 
f_cov = figure;
set(f_cov , 'Position', [680,314,1003,663]);
imagesc(Cov);
xticks(1:33);
yticks(1:33);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
xticklabels(para_label )
yticklabels(para_label)
title('Covariance matrix of parameters in \Theta_2 ','FontSize',15)
colorbar
colormap('turbo')
saveas(f_cov , fullfile(newdir  , ['FigS14_S3'  '.tif'] ))
PlotAllIntermediates_UpTodUMP(Para_trun, Para_mean_col )



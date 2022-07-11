function  f2 = MyStaAnalysis_dNTP_DSB(f2)
newdir_parent= 'Plot_Cellular_UpTodUMP_DSB';   % Your destination folder 
newdir  = fullfile(newdir_parent , 'DSB'  );
if ~isfolder( newdir   ) 
   mkdir( newdir   ) ;
else 
  rmdir( newdir ,'s'  ) ;
  mkdir( newdir ) ;
end
load("Para_dNTP_DSB_col_2.mat", 'Para_col')
Para_col_for_analysis = Para_col;
dbstop if error
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
percent = 0.2;%0.3
num_sample = size(Para_col_for_analysis,1);
num_para = length(para_label);
Para_trun = Para_col_for_analysis( floor(num_sample*percent)+1:end,1:num_para);
Para_mean_col = zeros(num_para,1);
T_col =  cell(num_para ,2);
T_onlyest =  cell(num_para ,2);
T_col(:,2) = num2cell(para_label );
T_onlyest(:,2) = num2cell(para_label );
num_trun = size(Para_trun,1);
factor = 25; 
for iFig = 1: ceil(num_para/factor)
    f = figure;
    set(0, 'CurrentFigure', f)
    set(f, 'Position', get(0, 'Screensize'));
    for i = 1: factor 
        subplot(5,5, i)
        idx = factor *(iFig-1)+i;
        dbstop if error
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
        CI  = quantile(x, [0.25, 0.75]);
        %CI = round(CI,);
        %bar_label = strcat( ' bar', '{', para_label(i),'}'   );
       %bar_label = strcat( '(',para_label(idx ),')','_{avg}' );
        temp = strcat( num2str( Para_mean) ,'(', num2str(CI(1)), ',', num2str(CI(2)),')');
        T_col{idx,1} = temp;
        %text( ax.XTick(1),  (ax.YTick(end)+ax.YTick(end-1))/2  ,temp,'tex'  )
        hold off
        if  idx == num_para
            break
        end
    end
    sgtitle( 'Histograms for parameters in \Theta_3','FontSize', 15)
    saveas( f, fullfile(newdir  , ['hist'  num2str(iFig) '.tif']));

end
save CI_dNTP_DSB.mat T_col  T_onlyest

num_sample_afterburnin = size(Para_trun,1);
Cov = zeros(num_para, num_para );
for i = 1:num_sample_afterburnin
    Cov  = Cov + (Para_trun(i,:) - Para_mean_col')' * (Para_trun(i,:) - Para_mean_col');
end
Cov = 1/(num_sample_afterburnin-1)*Cov; 
f_cov = figure;
set(f_cov , 'Position', [680,314,1003,663]);
imagesc(Cov);
xticks(1: num_para);
yticks(1:33);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
xticklabels(para_label )
yticklabels(para_label)
title('Covariance matrix of parameters in \Theta_3 ','FontSize',15)
colorbar
colormap('turbo')
saveas(f_cov , fullfile(newdir  , ['FigS15_S3'  '.tif'] ))


f2 = PlotIntermediates_dNTP_DSB(f2, Para_trun, Para_mean_col )
 
end
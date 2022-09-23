close all
clear
portion_coll = [0.7; 0.3;0.3;0.07;0.2;0;0.15;0;0.15 ];

f2 = figure;
set(f2, 'Position', get(0, 'Screensize'));
f_cov = figure;
set(f_cov , 'Position',  get(0, 'Screensize'));

newdir_parent  = 'Plot_TGI';
if ~isfolder(newdir_parent  )
    mkdir( newdir_parent   ) ;
end
for Model_Index = 1:9
    %     newdir = fullfile(newdir_parent, strcat("Model",num2str(Model_Index)) ) ;   % Your destination folder
%     if ~isfolder( newdir )
%         mkdir( newdir  ) ;
%     else
%         delete(  fullfile(newdir  ,'*' )   )
%         mkdir( newdir  ) ;
%     end
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
    
    factor = 9;
    T_col =  cell(1,num_para);
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
            CI  = quantile(x, [0.25, 0.75]);
            %CI = round(CI,3);
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
        sgtitle(['Histogram for parameters in model'   num2str(Model_Index), '  under uniform prior distribution'],'FontSize', 18)
       % saveas( f, fullfile(newdir_parent  , ['hist'  '_Model' num2str(Model_Index) '.tif']));
        exportgraphics(f,fullfile(newdir_parent  , ['hist'  '_Model' num2str(Model_Index) '.tif']), 'Resolution', 300 )

    end
    T = table(T_col);
    
    num_sample_afterburnin = size(Para_trun,1);
    Cov = zeros(num_para, num_para );
    for i = 1:num_sample_afterburnin
        Cov  = Cov + (Para_trun(i,:) - Para_mean_col')' * (Para_trun(i,:) - Para_mean_col');
    end
    Cov = 1/(num_sample_afterburnin-1)*Cov;
    set(0, 'CurrentFigure',     f_cov )
    subplot(3,3, Model_Index)
    imagesc(Cov);
    xticks(1:length(para_label));
    yticks(1:length(para_label));
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
    xticklabels(para_label )
    yticklabels(para_label)
    title( ['Model' num2str(Model_Index)],  'FontSize',15)
    colorbar
    colormap('parula')
    
    
    
    %     C_diff_cell =  kinetics_TumorVolume_Treated(Para_mean_col,Model_Index);
    %     lss = vecnorm([C_diff_cell{1}; C_diff_cell{2}],2);
    Plot_allTGIModel(f2,Para_trun,Para_mean_col,Model_Index);
end
set(0, 'CurrentFigure',     f_cov )
sgtitle( 'Covariance matrix of parameters in \Theta_4','FontSize',15)
exportgraphics(f2,fullfile( newdir_parent,'AllTGI.tif'), 'Resolution', 300 )
%saveas(f_cov , fullfile(newdir_parent , ['Cov_tumor'  '.tif'] ))
%saveas( f2,  fullfile( newdir_parent,'AllTGI.tif'));


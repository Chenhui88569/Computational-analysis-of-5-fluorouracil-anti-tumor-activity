function  PlotAllIntermediates_UpTodUMP(theta_col,theta)
%%
%1. Aim:Plot the selected intermediates in one figure
%2. Outcome: The kinetics of each concerned intermediates occupys one plot 
%3.  The functions that are responsible for determinding the time course of the intermediates are @Plot_kinetics_Cellular3to10
%4. The data sets of each essential intermediates are restored in the excel sheet '.\DataSet_all.xlsx"; errorbar in:DNA, RNA, TS
%5. The prediction intervals are plotted for each model variable
%6. The input argument theta means theta_mean
%7. The input argument theta_col means paramter popluation. The dimension
%is num_smaple * num_para
%% 
%addpath('/Users/chenhuima/MatlabWorkShop/PKPD_5FU/compartmental analysis')
dbstop if error
myfun_plot = @Plot_kinetics_Cellular3to10; 
Table_dir = 'DataSet_all.xlsx';
[T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir) ;
num_sample = 500; %size(theta_col,1);
idx_shuffle = randperm(size(theta_col,1),num_sample);
theta_col = theta_col(idx_shuffle,:);
%% interstitial 5FU, intracellular 5FU and 5FU anabolism

%% create kinetics based on the time span of DNA and RNA 
t_RNA = cell2mat( T_allsheets_output_cell(4)  );
time_span_max = min( t_RNA ): 0.5: max(t_RNA); % Make curve smooth ; t_RNA is zero
num_TimePoint = length(time_span_max);

Kinetics_Pop_inter = zeros(num_sample, num_TimePoint);
Kinetics_Pop_intra = zeros(num_sample, num_TimePoint);
Kinetics_Pop_ana = zeros(num_sample, num_TimePoint);
Kinetics_Pop_DNA = zeros(num_sample, num_TimePoint);
Kinetics_Pop_RNA = zeros(num_sample, num_TimePoint);
Kinetics_Pop_FreeTS = zeros(num_sample, num_TimePoint);
Kinetics_Pop_dUMP = zeros(num_sample, num_TimePoint);

dir_log = dir('*.log');
if ~isempty( [dir_log.name ])
    delete('*.log')
end
parfor i = 1: num_sample
    task  = getCurrentTask;
    taskid = task.ID;
    filename = ['logfile_' num2str(taskid) '.log'];
    outfile = fopen(  filename , 'a');
    theta_curr = theta_col(i,:);
    [~,Cfit_FNUC]  =  myfun_plot(theta_curr, time_span_max);
    %Cfit_col = Cfit_FNUC(:,[3,4,5])';
    Kinetics_Pop_inter(i,:) = Cfit_FNUC(:,3)';  % Cfit_col(1,:) ;
    Kinetics_Pop_intra(i,:) = Cfit_FNUC(:,4)';  %  Cfit_col(2,:) ;
    Kinetics_Pop_ana(i,:) = Cfit_FNUC(:,5)'; % Cfit_col(3,:) ;
    Kinetics_Pop_RNA(i,:) =  Cfit_FNUC(:,6)';
    Kinetics_Pop_DNA(i,:) =  Cfit_FNUC(:,7)';
    Kinetics_Pop_FreeTS(i,:)  =  Cfit_FNUC(:,end)';
    Kinetics_Pop_dUMP(i,:)  =  Cfit_FNUC(:,8)';
    if rem(num_sample,50) == 0
        fprintf(  outfile ,'sample , %d \r', i);
    end
end

%% plot preparation
label_add =  append ("(", ["A","B","C","D","E","F","G","I","H"] ,")");
f2 = figure;
set(f2, 'Position', get(0, 'Screensize'));
%% Plot under FNUC
t_FNUC = cell2mat( T_allsheets_output_cell(3)  );
FNUC_idx = find( time_span_max < max( t_FNUC) ) ;
FNUC_idx = FNUC_idx(end) + 1;
[Tv_FNUC_mean,Cfit_FNUC_mean]  =  myfun_plot(theta,time_span_max(1:FNUC_idx ));
Kinetics_Mean_inter = Cfit_FNUC_mean(:,3);
CI = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_inter(:,1:FNUC_idx),  Kinetics_Mean_inter, num_sample    );
%interstitial 5FU

subplot(3,3,1)
plot(Tv_FNUC_mean, Kinetics_Mean_inter ,'color',[0.85,0.33,0.10],'LineWidth',2); %mean kinetics 
hold on
[ha hb hc] = shadedplot( Tv_FNUC_mean', transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
box off
xlabel('Time(min)')
ylabel('amount(pmol/mg tissue)')
title(strcat( label_add(1), "  5FU in interstitial fluid"),'FontSize',15)

Kinetics_Mean_intra = Cfit_FNUC_mean(:,4);
CI = TimeCourses_ConfidenceInterval(...
   Kinetics_Pop_intra(:,1:FNUC_idx),  Kinetics_Mean_intra, num_sample    );
%intracellular 5FU

subplot(3,3,2)
c_Intra = cell2mat( C_allsheets_output_cell(2)  ); 
plot(t_FNUC,c_Intra,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]);%data
hold on
plot(Tv_FNUC_mean, Kinetics_Mean_intra,'LineWidth',2)
hold on
[ha hb hc] = shadedplot( Tv_FNUC_mean', transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
hold off; box off 
xlabel('Time(min)')
ylabel('amount(pmol/mg tissue)')
title( strcat( label_add(2), "  Intracellular 5FU") ,'FontSize',15 )

%anabolites
subplot(3,3,3)
Kinetics_Mean_ana = Cfit_FNUC_mean(:,5);
CI = TimeCourses_ConfidenceInterval(...
   Kinetics_Pop_ana(:,1:FNUC_idx),  Kinetics_Mean_ana, num_sample    );
c_FNUC = cell2mat( C_allsheets_output_cell(3)  ); 
p1 = plot(t_FNUC, c_FNUC,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]);%data
hold on
p2 = plot(Tv_FNUC_mean, Kinetics_Mean_ana, 'LineWidth',2);%model
[ha hb hc] = shadedplot( Tv_FNUC_mean', transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
hold off; box off
xlabel('Time(min)')
ylabel('amount(pmol/mg tissue)')
title(strcat( label_add(3), "  5FU anabolites") , 'FontSize',15 )


%% plot under RNA and DNA
t_RNA = cell2mat( T_allsheets_output_cell(4)  );
c_RNA = cell2mat( C_allsheets_output_cell(4)  );
[Tv_genome,Cfit_genome]  =  myfun_plot(theta, time_span_max); 
RNA_errorbar = [c_RNA(1:2);   [ 4.528;5.881;1.249]*10^3* 2*10^-3  ;  c_RNA(6) ;  [0.145 ;0.035 ]*10^3* 2*10^-3   ]; 
RNA_errorbar = abs(c_RNA- RNA_errorbar);

CI = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_RNA,  Cfit_genome(:,6), num_sample    );
%RNA
%figure
 subplot(3,3,4)
 e1 = errorbar( t_RNA,  c_RNA,RNA_errorbar,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data
 hold on 
 p2 = plot( t_RNA,  c_RNA,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410], 'MarkerEdgeColor',[0 0.4470 0.7410]  );%data
 hold on
 p1 = plot(Tv_genome, Cfit_genome(:,6),'color',[0.8500 0.3250 0.0980],'LineWidth',2);%model
 hold on
 [ha hb hc] = shadedplot(Tv_genome', transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
 hold off; box off
 xticks([24;48;72;168;240 ]*60)
 xticklabels({'1','2','3','7','10'} )
 xlabel('Time(day)')
 ylabel('amount (pmol/mg tissue)')
 title(strcat( label_add(4),"  F-RNA"), 'FontSize',15)
 legend([e1 p2 p1 ha],{'Data with errorbar','Data without errorbar', 'Model Simulation','','Prediction interval'},'FontSize',13,'Location','northeastoutside')
 legend boxoff 
 
 %DNA
 CI = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_DNA, Cfit_genome(:,7), num_sample    );
 subplot(3,3,5) 
 t_DNA = cell2mat( T_allsheets_output_cell(5)  );
 c_DNA = cell2mat( C_allsheets_output_cell(5)  ); 
 DNA_errorbar  = [ c_DNA(1:2);   16.626*0.2*10^-3 ; c_DNA(4) ; 13.349*0.2*10^-3;  c_DNA(6:end)      ];
 DNA_errorbar = abs(c_DNA- DNA_errorbar);
 data_plt= errorbar(t_DNA, c_DNA,DNA_errorbar, 'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data
 hold on
 hlp =  plot(Tv_genome, Cfit_genome(:,7), 'LineWidth',2 );%model
 hold on 
 [ha hb hc] = shadedplot(Tv_genome', transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
 hold off; box off 
 xlabel('Time(day)')
 ylabel('amount (pmol/mg tissue)')
 title( strcat( label_add(5), "  F-DNA") ,'FontSize',15)
 xticks([24;48;72;168;240]*60)
 xticklabels({'1','2','3','7','10'} )
 
 %% plot under TS inhibition
 c_complex = cell2mat( C_allsheets_output_cell(6)  ); %TS-FdUMP complex
 c_FreeTS = cell2mat( C_allsheets_output_cell(7)  );
 c_TotalTS = c_complex+ c_FreeTS ;
 c_FreeTS_errorbar = [11.833; 0; 0; 14.861; 24.984 ;34.700]*10^-3;
% c_FreeTS_errorbar = abs( c_FreeTS_errorbar- c_FreeTS);
 c_TotalTS_errorbar = [10.847 ;45.710;23.518;44.317;39.200;44.673]*10^-3;
 %c_TotalTS_errorbar = abs(c_TotalTS_errorbar -  c_TotalTS);
 c_FreeTSPercentage_errorbar=abs( c_FreeTS_errorbar./ c_TotalTS_errorbar -  c_FreeTS./ c_TotalTS  );
 t_TS = cell2mat( T_allsheets_output_cell(6)  );
 idx_TS = find(time_span_max < max(t_TS) );
 idx_TS = idx_TS (end) + 1 ;
 [Tv_TS,Cfit_TS]  =  myfun_plot( theta, time_span_max(1:idx_TS));
 [T_data_xtick_h, T_data_xlabel_h] = XtickLabelCreation( t_TS,t_TS/60 );
 
 

 %FreeTS
 CI = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_FreeTS(:,1: idx_TS ),Cfit_TS(:,end), num_sample    );
 %figure
 subplot(3,3,6) 
 errorbar( t_TS, c_FreeTS./c_TotalTS , c_FreeTSPercentage_errorbar,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data 
 hold on
 plot(Tv_TS, Cfit_TS(:,end), 'LineWidth',2, 'DisplayName','model');%model 
 hold on 
 [ha hb hc] = shadedplot(Tv_TS', transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval

 xticks(T_data_xtick_h)
 xticklabels(T_data_xlabel_h)
 xlabel('Time(h)')
 ylabel('amount (pmol/mg tissue)')
 title(strcat( label_add(6), "  % Free TS"),'FontSize',15)
 box off
 %legend({'Data','Simulation'},'FontSize',10) 


 %% plot under dUMP
 t_dUMP = cell2mat( T_allsheets_output_cell(8)  );
 c_dUMP = cell2mat( C_allsheets_output_cell(8)  );
 idx_dUMP = find(time_span_max < max(t_dUMP))+1;
 idx_dUMP =  idx_dUMP(end) + 1;
 [Tv_dUMP,Cfit_dUMP]  =  myfun_plot( theta, time_span_max(1: idx_dUMP));
 [T_data_xtick_h, T_data_xlabel_h] = XtickLabelCreation(  t_dUMP, t_dUMP/60 );
 
 CI = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_dUMP(:,1:idx_dUMP), Cfit_dUMP(:,8), num_sample    );
 
%figure 
subplot(3,3,7)
p1 = plot(  t_dUMP, c_dUMP  ,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data
hold on
p2 = plot(Tv_dUMP, Cfit_dUMP(:,8), 'LineWidth',2, 'DisplayName','model');%model
hold on
[ha hb hc] = shadedplot(Tv_dUMP', transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
hold off; box off
xticks(T_data_xtick_h)
xticklabels(T_data_xlabel_h)
xlabel('Time(h)')
ylabel('amount (pmol/mg tissue)')
title( strcat( label_add(7), "  dUMP") ,'FontSize',15)
 

f2 =  MyStaAnalysis_dNTP_DSB(f2);
newdir = 'Plot_Cellular_UpTodUMP_DSB';   % Your destination folder
saveas( f2, fullfile(newdir , 'AllCellular.tif'));

% %newdir  =strcat(FolderName , ['/', 'q', num2str(idx)  ]   );
% if ~isfolder( newdir ) 
%    mkdir( newdir  ) ;
% else 
%   delete( fullfile(newdir,'*' )  ) ;
% end
% FigList = findobj(allchild(0), 'flat', 'Type', '%figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   ohf = findobj( FigHandle); 
%   ohf_axes = findobj(ohf , 'Type','axes' );
%   if size(ohf_axes,1) == 2  % X and Y axes
%       FigHandle.CurrentAxes.TitleFontSizeMultiplier = 2;
%   end
%   %ohf = findobj(gcf)
% %   FigHandle = gca;
% %   TitleFontSizeMultiplier  % Scale factor for title font size
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   %set(FigHandle, 'Position', get(0, 'Screensize'));
end

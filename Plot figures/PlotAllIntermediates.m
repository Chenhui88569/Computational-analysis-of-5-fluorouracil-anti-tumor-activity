%1. Aim:Plot the selected intermediates in one figure
%2. Outcome: The kineictics of each concerned intermediates occupys one subplot 
%3.  The functions that are responsible for determinding the time course of
%  the intermediates are @kinetics_uptoFreeTS_SingleDose ,
%  @kinetics_dNTP_plot, @kinetics_DSB_plot
%4. The data sets of each essential intermediates are restored in the excel
% sheet '.\DataSet_all.xlsx';
% errorbar in:DNA, RNA, TS
clear
close all
myfun_plot = @kinetics_uptoFreeTS_SingleDose; %startfrom PK
Table_dir = 'DataSet_all.xlsx';
[T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir) ;
  
%% interstitial 5FU, intracellular 5FU and 5FU anabolism
t_FNUC = cell2mat( T_allsheets_output_cell(3)  );

[Tv_FNUC,Cfit_FNUC]  =  myfun_plot( [min(t_FNUC), max(t_FNUC)]);

figure 
%interstitial 5FU
%subplot(3,3,1)
plot(Tv_FNUC, Cfit_FNUC(:,3),'color',[0.85,0.33,0.10],'LineWidth',2)
box off
xlabel('Time(min)')
ylabel('amount(pmol/mg tissue)')
%title('5FU in interstitial fluid')

%intracellular 5FU
figure
%subplot(3,3,2)
c_Intra = cell2mat( C_allsheets_output_cell(2)  ); 
plot(t_FNUC,c_Intra,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]);%data
hold on
plot(Tv_FNUC, Cfit_FNUC(:,4),'LineWidth',2)
hold off; box off 
xlabel('Time(min)')
ylabel('amount(pmol/mg tissue)')
%title('Intracellular 5FU')

%anabolites
%subpsot(3,3,3)
figure
c_FNUC =cell2mat( C_allsheets_output_cell(3)  ); 
plot(t_FNUC, c_FNUC,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]);%data
hold on
plot(Tv_FNUC, Cfit_FNUC(:,5), 'LineWidth',2);%model
hold off; box off
xlabel('Time(min)')
ylabel('amount(pmol/mg tissue)')
%title('5FU anabolites')


%% RNA and DNA
 t_RNA = cell2mat( T_allsheets_output_cell(4)  );
 c_RNA = cell2mat( C_allsheets_output_cell(4)  );
[Tv_genome,Cfit_genome]  =  myfun_plot( [min( t_RNA ), max(t_RNA)]); 
RNA_errorbar = [c_RNA(1:2);   [ 4.528;5.881;1.249]*10^3* 2*10^-3  ;  c_RNA(6) ;  [0.145 ;0.035 ]*10^3* 2*10^-3   ]; 
RNA_errorbar = abs(c_RNA- RNA_errorbar);
%RNA
figure
 %subplot(3,3,4)
 errorbar( t_RNA,  c_RNA,RNA_errorbar,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data
 hold on
 plot(Tv_genome,Cfit_genome(:,6), 'LineWidth',2);%model
 hold off; box off
 xticks([24;48;72;168;240 ]*60)
 xticklabels({'1','2','3','7','10'} )
 xlabel('Time(day)')
 ylabel('amount (pmol/mg tissue)')
 %title('F-RNA')
 
 %DNA
 figure
 %subplot(3,3,5) 
 t_DNA = cell2mat( T_allsheets_output_cell(5)  );
 c_DNA = cell2mat( C_allsheets_output_cell(5)  ); 
 DNA_errorbar  = [ c_DNA(1:2);   16.626*0.2*10^-3 ; c_DNA(4) ; 13.349*0.2*10^-3;  c_DNA(6:end)      ];
 DNA_errorbar = abs(c_DNA- DNA_errorbar);
 data_plt= errorbar(t_DNA, c_DNA,DNA_errorbar, 'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data
 hold on
 hlp =  plot(Tv_genome,Cfit_genome(:,7), 'LineWidth',2 );%model
 hold off; box off 
 xlabel('Time(day)')
 ylabel('amount (pmol/mg tissue)')
 %title('F-DNA')
 xticks([24;48;72;168;240]*60)
 xticklabels({'1','2','3','7','10'} )
 
 %% TS inhibition
 c_complex = cell2mat( C_allsheets_output_cell(6)  ); %TS-FdUMP complex
 c_FreeTS = cell2mat( C_allsheets_output_cell(7)  );
 c_TotalTS = c_complex+ c_FreeTS ;
 c_FreeTS_errorbar = [11.833; 0; 0; 14.861; 24.984 ;34.700]*10^-3;
% c_FreeTS_errorbar = abs( c_FreeTS_errorbar- c_FreeTS);
 c_TotalTS_errorbar = [10.847 ;45.710;23.518;44.317;39.200;44.673]*10^-3;
 %c_TotalTS_errorbar = abs(c_TotalTS_errorbar -  c_TotalTS);
 c_FreeTSPercentage_errorbar=abs( c_FreeTS_errorbar./ c_TotalTS_errorbar -  c_FreeTS./ c_TotalTS  );
 t_TS = cell2mat( T_allsheets_output_cell(6)  );
 [Tv_TS,Cfit_TS]  =  myfun_plot(  [0, max( t_TS)]);
 [T_data_xtick_h, T_data_xlabel_h] = XtickLabelCreation( t_TS,t_TS/60 );

 %FreeTS
 figure
 %subplot(3,3,6) 
 errorbar( t_TS, c_FreeTS./c_TotalTS , c_FreeTSPercentage_errorbar,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data 
 hold on
 plot(Tv_TS, Cfit_TS(:,end), 'LineWidth',2, 'DisplayName','model');%model 
 xticks(T_data_xtick_h)
 xticklabels(T_data_xlabel_h)
 xlabel('Time(h)')
 ylabel('amount (pmol/mg tissue)')
 %title('% Free TS')
 box off
 %legend({'Data','Simulation'},'FontSize',10) 


 %% dUMP
 t_dUMP = cell2mat( T_allsheets_output_cell(8)  );
 c_dUMP = cell2mat( C_allsheets_output_cell(8)  );
 [Tv_dUMP,Cfit_dUMP]  =  myfun_plot(  [0, max(t_dUMP )]);
 [T_data_xtick_h, T_data_xlabel_h] = XtickLabelCreation(  t_dUMP, t_dUMP/60 );
 
figure 
%subplot(3,3,7)
plot(  t_dUMP, c_dUMP  ,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410]  );%data 
 hold on
 plot(Tv_dUMP, Cfit_dUMP(:,8), 'LineWidth',2, 'DisplayName','model');%model
 hold off; box off
 xticks(T_data_xtick_h)
 xticklabels(T_data_xlabel_h)
 xlabel('Time(h)')
 ylabel('amount (pmol/mg tissue)')
 %title('dUMP')
 
 
 %% dNTP pool imbalance
 figure
 %subplot(3,3,8)
Anabolites_UnitConversion =10^6/130.077  ;
day2min = 60*24;
Dose_info = [7*day2min	 0	7*day2min	 92.9936306*Anabolites_UnitConversion]; % A week
t_dNTP =cell2mat( T_allsheets_output_cell(9)  ); 
t_dNTP_hour = t_dNTP/60; addpoint = [2,4,8,9];
[t_dNTP_xtick, t_dNTP_xlabel] = XtickLabelCreation(t_dNTP,t_dNTP_hour );
           
theta_dNTP = [0.21351;32.50090;3.71407;73.34474;3.77027;0.30262;3.33454;1.35311;0.27614;35.66968;6.04530;0.07853];
myfun_dNTP_plot = @kinetics_dNTP_plot;
Dose_info(1) = t_dNTP(end);
[Tv,Cfit]  =  myfun_dNTP_plot(theta_dNTP,  Dose_info);  

Current_dNTP_M = 'alldNTP';
C_dNTP_measure_pertubation = dNTP_data_processing(Current_dNTP_M);
C_data = C_dNTP_measure_pertubation(:,1);
plot( t_dNTP(~ismember(t_dNTP_hour,addpoint) ), C_data(~ismember(t_dNTP_hour,addpoint ) ),'o', 'MarkerSize',6.5,'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'data' );%data
hold on
hlp =  plot(Tv, Cfit(:,1), 'LineWidth',2 );%model
hold off; box off 
xticks(t_dNTP_xtick)
xticklabels(t_dNTP_xlabel )
xlabel('Time(h)')
ylabel('pertubation')
%title('dNTP pool imbalance')

%%  DSB
t_DSB =cell2mat( T_allsheets_output_cell(10)  ); 
c_DSB = cell2mat( C_allsheets_output_cell(10)  ); 
t_DSB_hour = t_DSB/60 ;
[t_DSB_xtick, t_DSB_xlabel] = XtickLabelCreation(t_DSB,t_DSB_hour );
Dose_info(1) =t_DSB(end);  
theta_DSB = [ 2577.17019 10.50768 151.94553 10.06779 194.02545 0.12550 2.79398];
[Tv,Cfit]  =   kinetics_DSB_plot(theta_DSB ,  Dose_info );
figure
%subplot(3,3,9)              
plot(t_DSB(~ismember(t_DSB_hour, [32  64]) ), c_DSB(~ismember(t_DSB_hour, [32  64]))  ,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'data');%data
hold on
plot(Tv, Cfit, 'LineWidth',2  );%model 
hold off; box off
xticks(t_DSB_xtick)
xticklabels(t_DSB_xlabel )
xlabel('Time(h)')
ylabel('γ-H2AX foci intensity(thousands)')
%title('Double strand break generation')
              
 
idx = 1;
newdir = '../plot/Intermediates_plot';   % Your destination folder
%newdir  =strcat(FolderName , ['/', 'q', num2str(idx)  ]   );
if ~isfolder( newdir ) 
   mkdir( newdir  ) ;
else 
  delete( fullfile(newdir,'*' )  ) ;
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  saveas( FigHandle, fullfile(newdir , [FigName '.png']));
end

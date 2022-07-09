clear
close all
%addpath('../DataSet_TGI')
dbstop if error
warning('off')
%Parameters fall into five category
%The one element means the parameters in certain segment is involved
Anabolites_flag_col = [0,0,1,0,1,0,0,0,0];   Anabolites_flag_col([5,7]) =[];
TSinhibion_flag_col = [1,1,1,1,0,0,1,0,1];   TSinhibion_flag_col([5,7]) = [];
dNTP_flag_col =       [0,0,0,0,0,1,0,1,0]; dNTP_flag_col([5,7]) = [];
DSB_flag_col =        [0,0,0,1,1,1,1,0,0];   DSB_flag_col([5,7]) = [];
TGI_flag_Col =        [0,0,0,0,0,0,0,1,1];   TGI_flag_Col([5,7]) = [];
TSInhibition_TestedPara_col = {["alpha", "59","95","09"],["G_0","cat","08",  "59","95","09"], [ "59","95","09"],...
   ["59","95","09"], [],[],["59","95","09"],[],["59","95","09"]};
TSInhibition_TestedPara_col([5,7]) = [];
Anabolites_TestedPara_col = {[],[],["V_max_54","K_m_54"],[],["V_max_54","K_m_54"],[],[],[],[]};
Anabolites_TestedPara_col([5,7]) = [];
DSB_TestedPara_col = {[],[],[],"lag_DSB","lag_DSB","lag_DSB",["V_HR","K_HR"],[],[]};
DSB_TestedPara_col([5,7]) = [];
TGI_TestedPara_col = {[],[],[],[],[],[],[],"lag_tumor_TGI",["E_max_damage","EC_50_damage"]};
TGI_TestedPara_col([5,7]) = [];
%A = figure('Position',get(0, 'Screensize'));

newdir = 'plot_Resistance' ;   % Your destination folder
if ~isfolder( newdir ) 
   mkdir( newdir  ) ;
else
   rmdir(newdir,'s');
   mkdir( newdir  ) ;
end

filename_ParaChange = ['logfile_'  'ParaChange' '.log'];
delete(filename_ParaChange)
%% The final product I want to get is that a figure has 7 panels. Each subpanel has 3 panels.
f_resis = figure;
set(f_resis,'Position' ,get(0,'Screensize'));

num_cases = length(Anabolites_flag_col); % make sure it's 7
parpool(num_cases);
for i = 1:num_cases 
    %count = i-3*(batch_num-1);  %for plot,starting from 1 for each batch
    TSInhibition_flag = TSinhibion_flag_col(i);
    Anabolites_flag = Anabolites_flag_col(i);
    dNTP_flag = dNTP_flag_col(i);
    DSB_flag = DSB_flag_col(i); 
    TGI_flag = TGI_flag_Col(i); 
    TSInhibition_TestedPara  = TSInhibition_TestedPara_col(i);
    TSInhibition_TestedPara  = [TSInhibition_TestedPara{:}];
    Anabolites_TestedPara = Anabolites_TestedPara_col(i);
    Anabolites_TestedPara  = [Anabolites_TestedPara{:}]; 
    DSB_TestedPara  = DSB_TestedPara_col{i};
    TGI_TestedPara  = TGI_TestedPara_col{i};
    ResistanceAnalysis(f_resis ,newdir,i, Anabolites_flag, Anabolites_TestedPara,...
        TSInhibition_flag ,TSInhibition_TestedPara,...
        dNTP_flag,DSB_flag,  TGI_flag,DSB_TestedPara, TGI_TestedPara )
end
delete(gcp('nocreate'))
saveas( f_resis, fullfile(newdir , 'AllResis.tif'));
saveas( f_resis, fullfile(newdir , 'AllResis.fig'));

% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   saveas( FigHandle, fullfile(newdir , ['Case_', FigName '.png']));
% end


function ResistanceAnalysis(f_resis, newdir,Casenum, Anabolites_flag, Anabolites_TestedPara,...
    TSInhibition_flag ,TSInhibition_TestedPara,...
    dNTP_flag,DSB_flag, TGI_flag,DSB_TestedPara, TGI_TestedPara ) 

UnitConversion = [  10^6/130.077, 24*60,10^-3 ,7*24*60   ];
UnitConversion_cell = num2cell(UnitConversion  );
[Plasma_UnitConversion, day2min  ,mm3tocm3  ,week2min  ] = UnitConversion_cell{:};
Table_dir = 'DataSet_all.xlsx';
[T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir);
lag_tumorNTP = cell2mat(T_allsheets_output_cell(9) ) ;
t_TS = cell2mat(T_allsheets_output_cell(6) ) ;
t_DSB   = cell2mat(T_allsheets_output_cell(10) ) ; %last data point
t_anabolites = t_TS ;
Data_array_Tv = importdata('DataSet_TGI/50mgkg_EveryTheOtherDay_LS174T_resistant.mat');
Dose_info = [15*60*24  5  48*60  93/2];
theta_TGI_50_sensitive = [1339.31457867841;104.059993350524;0.841096217568630;73.9835830977462;0.000227131913701615;11.2628519108548;4.57756840800819e-06];
theta_TGI_50_resistance = [1300.244592
72.20112535
1.524956809
49.65603054
0.0002737007422
7.564421928
5.51e-06];

Tv_initial = 125*mm3tocm3;
t_Tumor = Data_array_Tv(:,1); %The prediction doesn't proceed beyond the experimental data
%% assign parameter value to the submodels
load('Para.mat','Para_Info');
Para_value = Para_Info.Para_Value;
theta_cellular = Para_value(1:32);
%theta_FNUC = [25.17731647	2321.536075];
V_max_54 = theta_cellular(12);
K_m_54 = theta_cellular(13);
theta_TS = [theta_cellular(26:end);2.0212;1.0345;0.0186]; 
theta_TS_cell = num2cell(theta_TS);
[k_95, k_59 ,k_09, K_dUMP, G_0, k_cat, k_08, alpha, k_d, TS_0] = theta_TS_cell{:}; 

theta_dNTP = [ Para_value(33:33+11); 0.15 ];
theta_dNTP_cell = num2cell(theta_dNTP);
[k1, k2,k3,k4,k5,k6,k7,k8,k9, k_B,k_A,k10 ,gamma_dNTP ] = theta_dNTP_cell{:};

theta_DSB = Para_value(33+11+1:33+18);
theta_DSB_cell = num2cell(theta_DSB);
[lag_DSB , V_dNTP, K_dNTP, V_HRR, K_HRR, k_i,k_0] =  theta_DSB_cell{:}; 

theta_TGI_cell = num2cell(theta_TGI_50_sensitive);
[lag_tumor, IC50,E_max_damage, EC_50 ,lambda_g , P_max, lambda_d ] = theta_TGI_cell{:}; 

Concerned_Para_original= [ V_max_54,K_m_54,G_0, k_cat, k_08, k_95,k_59 ,k_09, k5,lag_DSB,alpha, V_HRR, K_HRR, lag_tumor,E_max_damage, EC_50 ] ;
Elongation = 10*60;%additional 10 hours
DSB_TimeSpan  = min(t_DSB  ) : 20: max( t_DSB  )+Elongation;
Dose_info_DSB = Dose_info; 
if DSB_flag==1
    %When the effects of lag_DSB are examined, the end of time span/final
    %integration time should be elongated to display the full time profile. However, the last part of requested elongated
    %time would eventually lie outside the domain of treatment course.
    t_DSB(end) = t_DSB(end) + 20 *day2min;
    DSB_TimeSpan = min(t_DSB  ) : 100: max( t_DSB  );
    Dose_info_DSB(1)  = max( t_DSB  );
end
%the duration of the treatment should be larger than the simulaton time
%for the species 
Dose_info_Tv= Dose_info;
Dose_info_Tv(1) = 2*week2min+2*day2min; %make sure the that time span does not surpass the point near the equilibrium point no more than 3 weeks
TumorVolume_TimeSpan  = min( t_Tumor) :200:Dose_info_Tv(1)  ;
TumorVolume_Timelength =  size(TumorVolume_TimeSpan,2  );
DSB_Timelength =  size( DSB_TimeSpan,2  );
dNTP_TimeSpan  = min(lag_tumorNTP):20: max(lag_tumorNTP)+Elongation;
dNTP_Timelength = size(dNTP_TimeSpan,2  );
TSInhibition_TimeSpan  = min(t_TS):20: max(t_TS)+Elongation;
TSInhibition_Timelength = size(TSInhibition_TimeSpan,2  );
Anabolites_TimeSpan  = min(t_anabolites ):20: max(t_anabolites ) ;
Anabolites_Timelength = size(Anabolites_TimeSpan ,2  );

%V_max_54,K_m_54,G_0, k_cat,k08,k_95,k_59 ,k_09, k5,lag_DSB  


filename_ParaChange = ['logfile_'  'ParaChange' '.log'];
outfile_ParaChange = fopen(  filename_ParaChange , 'a');



Casenum_load = Casenum ;

outdir = ['outputfile_case' num2str(Casenum_load) '/'];
if isfolder(outdir)
   load([outdir  'Population_', num2str(Casenum_load),'.mat'], 'ParaPopulation_permutated','TumorVolume_PopulationKinetics','Anabolites_PopulationKinetics','DSB_PopulationKinetics');
else
    load(['Population_', num2str(Casenum_load),'.mat'], 'ParaPopulation_permutated','TumorVolume_PopulationKinetics','Anabolites_PopulationKinetics','DSB_PopulationKinetics');
end

TumorVolume_PopulationKinetics(TumorVolume_PopulationKinetics<0) = 0;
ParaPopulation_size = size(TumorVolume_PopulationKinetics,1);

[ Tsim_all,Ysim_all]  =  Resistance_kinetics_Tv( theta_cellular, theta_dNTP, theta_DSB, theta_TGI_50_sensitive,Concerned_Para_original,Tv_initial, Dose_info_DSB); 

% [T_UPtoFreeTS,Cv_UPToFreeTS] = Resistance_kinetics_uptoFreeTS(theta_cellular, Concerned_Para_original,Dose_info);
% [lag_tumorNTP,Cv_dNTP] = Resistance_kinetics_dNTP_plot(theta_cellular, theta_dNTP,Concerned_Para_original,Dose_info);
% [t_DSB,Cv_DSB] =   Resistance_kinetics_DSB_plot(theta_cellular, theta_dNTP,theta_DSB,Concerned_Para_original,Dose_info);

T_UPToFreeTS = Tsim_all{1};Tv_dNTP = Tsim_all{3} ; Tv_DSB = Tsim_all{4}; T_Tumor = Tsim_all{5};
Cv_ana = Ysim_all{1};  Cv_TS = Ysim_all{2};Cv_dNTP = Ysim_all{3} ; Cv_DSB = Ysim_all{4}; Cv_tumor = Ysim_all{5};
Cv_tumor = interp1(T_Tumor   , Cv_tumor, TumorVolume_TimeSpan     ,  'PCHIP'); %percentage free TS
Anabolites_original  = interp1( T_UPToFreeTS  ,Cv_ana , Anabolites_TimeSpan    ,  'PCHIP'); %percentage free TS
TSInhibition_original  = interp1( T_UPToFreeTS  ,Cv_TS , TSInhibition_TimeSpan    ,  'PCHIP'); %percentage free TS
dNTP_original =  interp1( Tv_dNTP ,Cv_dNTP(:,1) ,dNTP_TimeSpan  ,  'PCHIP');
DSB_original =  interp1( Tv_DSB,Cv_DSB ,DSB_TimeSpan  ,  'PCHIP');

Concerned_Para_original_round = round(Concerned_Para_original,3);
ParaLabel =[ "V_{max,54}","K_{m,54}", "G_{0,dUMP}","k_{cat,dUMP}","k_{08,dUMP}", "k_{95,complex}","k_{59,complex}","k_{09,complex}",...
    "k_{5,dNTP}","lag_{DSB} ","alpha", "V_{max,HR}" ,"K_{max,HR" , "T_{d,Tv}","E_{max,damage}","EC_{50, damage}" ];%string array

%  Turn nummeric array to string array using string() instead of num2str()!!!!! Using num2str() will convert the entire numeric array into one character.

%----------------------------------
%  The following codes are used to summarize the analyze the simulation results of all the parameter samples. 
%
% `idx` points to the parameter set that can give rise to the simulation result closest to the resistant response. 
%----------------------------------

[CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN ] = Resistance_InferentialAnalysis(...
    TumorVolume_PopulationKinetics, ParaPopulation_size, ParaPopulation_permutated    );%CI: confidence interval
[Tsim_50_resistance,Ysim_50_resistance] =  Resistance_kinetics_Tv( theta_cellular, theta_dNTP, theta_DSB, theta_TGI_50_resistance,Concerned_Para_original, Tv_initial, Dose_info_Tv);
T_tumor_50_resistance = Tsim_50_resistance{end};
Cv_tumor_50_resistance = Ysim_50_resistance{end};
T_tumor_50_sensitive = TumorVolume_TimeSpan;
Cv_tumor_50_sensitive = Cv_tumor ;
[idx , ClosetTimeCourse] = CompareTimeCourse(T_tumor_50_resistance,Cv_tumor_50_resistance,...
    TumorVolume_PopulationKinetics ,TumorVolume_TimeSpan );
para_ClosetTimeCourse =  ParaPopulation_permutated(idx,:);

[T_tumor_close_all,Cv_tumor_50_close_all] =  Resistance_kinetics_Tv( theta_cellular, theta_dNTP, theta_DSB, theta_TGI_50_sensitive, para_ClosetTimeCourse, Tv_initial, Dose_info_Tv);
 T_tumor_close = T_tumor_close_all{end};
 Cv_tumor_50_close = Cv_tumor_50_close_all{end};
%fprintf("Case: %d\n", count+3*(batch_num-1))
%subplot(Case_num,3,3*(count-1)+1) 


%%
letters = append("(",["A","B", "C","D", "E","F","G"], ")");
%% plotting
set(0,'CurrentFigure', f_resis);
ax = subplot(7,3, (Casenum-1)*3+1);

%f = figure('Position',[1,399,1440,398]);
%subplot(1,3,1) 
if C_MAX == C_MIN
    plot( T_tumor_50_sensitive, Cv_tumor_50_sensitive  , 'LineWidth',2 ,'DisplayName' ,'original'    );
    legend
else
    f1= plot( T_tumor_50_resistance,  Cv_tumor_50_resistance ,'r', 'LineWidth',2    );
    hold on
    f2 = plot( T_tumor_50_sensitive,Cv_tumor_50_sensitive ,'b', 'LineWidth',2    );
    hold on 
    f3  = plot(  TumorVolume_TimeSpan, C_mean ,'--','Color',[0.93,0.69,0.13], 'LineWidth',2    );
    hold on 
    f4 = plot( T_tumor_close,Cv_tumor_50_close ,':','Color',[0.49,0.18,0.56], 'LineWidth',2    );
    hold on
    [ha hb hc] = shadedplot( TumorVolume_TimeSpan', transpose(CI(:,1)),transpose(CI(:,2)), 'c','none');
    %ha is a 1*2 area array ,  corresponds to two descriptive labels
    hold on
    max_legend   = CompareWithOriginalValue( para_MAX, Concerned_Para_original,ParaLabel); %return a cell array of character vector
    min_legend  = CompareWithOriginalValue( para_MIN, Concerned_Para_original,ParaLabel);
    Closest_legend = CompareWithOriginalValue( para_ClosetTimeCourse' , Concerned_Para_original,ParaLabel);
end
title([ letters(Casenum) ' Case ' num2str(Casenum) ' :Tumor volume'],'FontSize', 13')
ax.TitleHorizontalAlignment = 'left';
grid off ;box off
xlabel('Time(week)')
ylabel('fold change')
xticks( 0:week2min: fix( Dose_info(1) /week2min)*week2min )
xticklabels( (  num2cell( 0:1: fix( Dose_info(1) /week2min) )) )
% title("Tumor Volume kinetics"+newline + "original parameters: "+ OriginPara_Display )  
fprintf(outfile_ParaChange,"%d:Tumor Volume\n", Casenum );
fprintf(outfile_ParaChange, "\t Max: %s \n",  strjoin( max_legend) );
fprintf(outfile_ParaChange, "\t Min: %s \n",  strjoin( min_legend) );
fprintf(outfile_ParaChange, "\t Closest: %s \n",  strjoin( Closest_legend) );
% legend([f1,f2,f3,f4,ha] , {"Resistant tumor response" ,"Sensitive tumor response",...
%         'Mean time course',...
%         'Time course with \Theta_{closest}',...
%         ' ', '90% confidence interval'},...
%         'Location','northwestoutside','FontSize', 13)

  
% [T_UPtoFreeTS_closeTv,Cv_UPToFreeTS_closeTv] = Resistance_kinetics_uptoFreeTS(theta_cellular, para_ClosetTimeCourse,Dose_info);
% [t_DSB_closeTv,Cv_DSB_closeTv] =   Resistance_kinetics_DSB_plot(theta_cellular, theta_dNTP,theta_DSB,para_ClosetTimeCourse,Dose_info);


T_UPtoFreeTS_closeTv =  T_tumor_close_all{1}  ;Anabolites_closeTv = Cv_tumor_50_close_all{1};
Tv_DSB_closeTv = T_tumor_close_all{end-1}; Cv_DSB_closeTv  =  Cv_tumor_50_close_all{end-1};
Anabolites_closeTv  = interp1(T_UPtoFreeTS_closeTv,Anabolites_closeTv, Anabolites_TimeSpan    ,  'PCHIP'); %percentage free TS
DSB_closeTv =  interp1(  Tv_DSB_closeTv,Cv_DSB_closeTv ,DSB_TimeSpan  ,  'PCHIP');
%para_MAX and para_MIN are column vectors
[CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN ] = Resistance_InferentialAnalysis(...
    Anabolites_PopulationKinetics, ParaPopulation_size,  ParaPopulation_permutated     );%CI: confidence interval
[lag_tumorata_xtick, lag_tumorata_xlabel] = XtickLabelCreation( t_anabolites  , t_anabolites /60 ); %lag_tumorata , lag_tumorata/24/60

subplot(7,3, (Casenum-1)*3+2)
if C_MAX == C_MIN
    plot( Anabolites_TimeSpan, Anabolites_original   , 'LineWidth',2,'DisplayName' ,'original'   );
else
    f1 = plot(  Anabolites_TimeSpan, Anabolites_original,'b', 'LineWidth',2    );
    hold on 
    f = plot( Anabolites_TimeSpan, C_mean,'--' ,'color',[0.93,0.69,0.13],'LineWidth',2    );
    hold on
    f2 =  plot( Anabolites_TimeSpan, Anabolites_closeTv, ':','color',[0.49,0.18,0.56], 'LineWidth',2    );
    hold on
    [ha hb hc] = shadedplot( Anabolites_TimeSpan', transpose(CI(:,1)),transpose(CI(:,2)), 'c','none');
    hold off
    max_legend   = CompareWithOriginalValue( para_MAX, Concerned_Para_original,ParaLabel);
    min_legend  = CompareWithOriginalValue( para_MIN, Concerned_Para_original,ParaLabel);
end
% legend([f1,f,f2,ha] , {"Simiulated anabolites in sensitive cell line" ,...
%         'Mean time course',...
%         'Time course with parameters closest to resistant cell line',...
%         ' ', '90% confidence interval'},...
%         'Location','northeast','FontSize', 10)
title(['Case ' num2str(Casenum) ,' :Anabolites'],'FontSize', 13')
grid off ;box off
xlabel('Time(h)')
ylabel('pmol/mg')
xticks( lag_tumorata_xtick )
xticklabels(lag_tumorata_xlabel )
fprintf(outfile_ParaChange,"%d: Anabolites\n" ,Casenum);
fprintf(outfile_ParaChange,"\t Max: %s \n",  strjoin( max_legend) );
fprintf(outfile_ParaChange,"\t Min: %s \n",  strjoin( min_legend) );
fprintf(outfile_ParaChange,"\t Closest: %s \n",  strjoin( Closest_legend) );


subplot(7,3, (Casenum-1)*3+3)
%subplot(1,3,3)
[lag_tumorata_xtick, lag_tumorata_xlabel] = XtickLabelCreation(t_DSB ,t_DSB /60 ); %lag_tumorata , lag_tumorata/24/60
[CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN ] = Resistance_InferentialAnalysis(...
    DSB_PopulationKinetics, ParaPopulation_size,   ParaPopulation_permutated     );%CI: confidence interval
C_mean(C_mean<0) = 0;
DSB_closeTv(DSB_closeTv<0) = 0;
len_mean = length(C_mean) ;
if C_MAX == C_MIN
    plot(  DSB_TimeSpan, DSB_original  , 'LineWidth',2,'DisplayName' ,'original'   );
else
    f1 = plot( DSB_TimeSpan, DSB_original ,'b', 'LineWidth',2    );
    hold on
    f = plot(  DSB_TimeSpan(1:len_mean),C_mean,'--','color',[0.93,0.69,0.13], 'LineWidth',2    ); 
    hold on
    f2 =  plot( DSB_TimeSpan, DSB_closeTv,':' ,'color',[0.49,0.18,0.56], 'LineWidth',2    ); 
    [ha hb hc] = shadedplot( DSB_TimeSpan' ,transpose(CI(:,1)),transpose(CI(:,2)), 'c','none');
    hold off
    max_legend   = CompareWithOriginalValue( para_MAX, Concerned_Para_original,ParaLabel);
    min_legend  = CompareWithOriginalValue( para_MIN, Concerned_Para_original,ParaLabel);
end
grid off ;box off

% legend([f1,f,f2,ha] , {"Simiulated DSB  in sensitive cell line" ,...
%         'Mean time course',...
%         'Time course with parameters closest to resistant cell line',...
%         ' ', '90% confidence interval'},...
%         'Location','northwest','FontSize', 10)

if DSB_flag ==1
    xlabel('Time(week)')
    ylabel('count(thousands)')
    xticks( 0:week2min: fix( Dose_info_DSB(1) /week2min)*week2min )
    xticklabels( (  num2cell( 0:1: fix( Dose_info_DSB(1) /week2min) )) )
else
    xlabel('Time(day)')
    ylabel('count(thousands)')
    xticks( min(t_DSB) :24*60: max(t_DSB) )
    xticklabels( num2cell(  min(t_DSB) : 1  : max(t_DSB)/60/24 )) 
end
fprintf(outfile_ParaChange,"%d: DSB\n", Casenum );
fprintf(outfile_ParaChange, "\t Max: %s \n",  strjoin( max_legend) );
fprintf(outfile_ParaChange, "\t Min: %s \n",  strjoin( min_legend) );
fprintf(outfile_ParaChange, "\t Closest: %s \n",  strjoin( Closest_legend) );
title(['Case ' num2str(Casenum) ,' :DSB'],'FontSize', 13')

%saveas( f, fullfile(newdir , ['Case_', num2str(Casenum) '.tif']));
end



 

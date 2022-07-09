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

ParaPopulation_size = 500; %total number of parameters that are sampled out.
%A = figure('Position',get(0, 'Screensize'));
Case_num = 1;%number of cases per batch
%batch_num = 2;

%rmdir('outputfile_case*')
for i = 1:9
    outdir = ['outputfile_case' num2str(i) '/'];
    if ~isfolder(outdir)
        mkdir(outdir);
    else 
        rmdir(outdir,'s');
        mkdir(outdir);
    end
    %dir_log = dir([outdir '*.log']);
    %if ~isempty( [dir_log.name])
      %  delete([outdir '*.log'])
   % end
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
    ResistanceAnalysis(outdir,i, Anabolites_flag, Anabolites_TestedPara,...
        TSInhibition_flag ,TSInhibition_TestedPara,...
        dNTP_flag,DSB_flag,  TGI_flag,DSB_TestedPara, TGI_TestedPara, ParaPopulation_size )
    delete(gcp('nocreate'))
end
%fullpath = strcat(folderpath, "\Allcase.png" );
%saveas(A,fullpath )
% newdir = fullfile('..',  'plot_Resistance' ) ;   % Your destination folder
% %newdir  =strcat(FolderName , ['/', 'q', num2str(idx)  ]   );
% if ~isfolder( newdir ) 
%    mkdir( newdir  ) ;
% else
%    delete(  fullfile(newdir  ,'*' )   )
% end
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   %saveas( FigHandle, fullfile(newdir , ['Case_', FigName '.png']));
% end


function ResistanceAnalysis(outdir,Casenum, Anabolites_flag, Anabolites_TestedPara,...
    TSInhibition_flag ,TSInhibition_TestedPara,...
    dNTP_flag,DSB_flag, TGI_flag,DSB_TestedPara, TGI_TestedPara, ParaPopulation_size ) 

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
theta_TGI_50_sensitive = [1115.97640747442
87.8008399139146
0.491458377564318
62.8294048603437
0.000231228259995630
8.96996460249570
4.13189231322516e-06
0.0115216673006455
0.0485615994786281];
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
if DSB_flag==1
    %When the effects of lag_DSB are examined, the end of time span/final
    %integration time should be elongated to display the full time profile. However, the requested elongated
    %time would lie outside the domain of treatment course, which
    %could be tackled either  by extraploation or by protracting the
    %duration of the treatment.
    t_DSB(end) = t_DSB(end) + 20 *day2min;
    DSB_TimeSpan = min(t_DSB  ) : 100: max( t_DSB  );
    Dose_info(1)  = max( t_DSB  );
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


%The time course variations have the same dimension,faciliating the further analysis.
Anabolites_PopulationKinetics = zeros(ParaPopulation_size , Anabolites_Timelength);
TSInhibtion_PopulationKinetics = zeros(ParaPopulation_size , TSInhibition_Timelength);
dNTP_PopulationKinetics = zeros(ParaPopulation_size , dNTP_Timelength );
DSB_PopulationKinetics = zeros(ParaPopulation_size , DSB_Timelength);
TumorVolume_PopulationKinetics =  zeros(ParaPopulation_size, TumorVolume_Timelength  );




%V_max_54,K_m_54,G_0, k_cat,k08,k_95,k_59 ,k_09, k5,lag_DSB  

%ParaPopulation_permutated = ParaPopulation_unpermutated;
% for i = 1: size(  Concerned_Para_original ,2) %shuffle the parameter population matrix
%     premuted_column = randperm( ParaPopulation_size  );
%     ParaPopulation_permutated(:,i) =  ParaPopulation_unpermutated( premuted_column ,i   );
% end


ParaPopulation_permutated = zeros(ParaPopulation_size, size(  Concerned_Para_original ,2));

parpool(180);
parfor i = 1:ParaPopulation_size %all the simulation should be have the same time span
    %disp(num2str(i));
    flag = 1;
    task  = getCurrentTask;
    taskid = task.Id;
    filename = [ outdir 'logfile_' num2str(taskid) '.log'];
    outfile = fopen(  filename , 'a');
    fprintf(   outfile, '%d th Case: Analysis on %d th set \n' ,  Casenum, i);
   % theta_variation = ParaPopulation_permutated;
%     [T_UPToFreeTS,Cv_UPToFreeTS] = Resistance_kinetics_uptoFreeTS(theta_cellular, theta_variation,Dose_info);
%     [Tv_dNTP,Cv_dNTP] = Resistance_kinetics_dNTP_plot(theta_cellular, theta_dNTP,theta_variation,Dose_info);
%     [Tv_DSB,Cv_DSB] =   Resistance_kinetics_DSB_plot(theta_cellular, theta_dNTP,theta_DSB,theta_variation,Dose_info);
%     [T_Tumor,Cv_tumor] =  Resistance_kinetics_Tv( theta_cellular, theta_dNTP, theta_DSB, theta_TGI_50_sensitive, theta_variation,Tv_initial, Dose_info); 
    count = 0;
    while flag == 1
        lastwarn('') ; %clear last warning.
        count =  count + 1;
        Para_permutated = CreateParaPopulation(Anabolites_flag, Anabolites_TestedPara,...
            TSInhibition_flag ,TSInhibition_TestedPara,...
            dNTP_flag,DSB_flag, TGI_flag,DSB_TestedPara, TGI_TestedPara, 1 );
        theta_variation = Para_permutated;
        [ Tsim_all,Ysim_all]  =  Resistance_kinetics_Tv( theta_cellular, theta_dNTP, theta_DSB, theta_TGI_50_sensitive, theta_variation,Tv_initial, Dose_info);
        T_UPToFreeTS = Tsim_all{1};Tv_dNTP = Tsim_all{3} ; Tv_DSB = Tsim_all{4};T_Tumor = Tsim_all{5};
        Cv_ana = Ysim_all{1};  Cv_TS = Ysim_all{2};Cv_dNTP = Ysim_all{3} ; Cv_DSB = Ysim_all{4}; Cv_tumor = Ysim_all{5};
        difference_Ysim = diff(Cv_tumor );
        sec_difference_Ysim = diff(difference_Ysim);
        if max(abs( difference_Ysim )) >= 10 || max(abs( sec_difference_Ysim )) >= 5 || max(abs(Cv_DSB))>150    
        %if ~isempty(lastwarn)
            flag = 1;
            fprintf(   outfile, '  %d th try failed \n' ,  count);
        else
            flag = 0;
            fprintf(   outfile, '  %d th try succeds \n' ,  count);           
        end
    end
    TumorVolume_PopulationKinetics(i,:) = interp1(T_Tumor , Cv_tumor   ,TumorVolume_TimeSpan,  'PCHIP');
    Anabolites_PopulationKinetics(i,:) =  interp1( T_UPToFreeTS ,  Cv_ana ,Anabolites_TimeSpan  ,  'PCHIP');
    TSInhibtion_PopulationKinetics(i,:) =  interp1( T_UPToFreeTS , Cv_TS ,TSInhibition_TimeSpan  ,  'PCHIP');
    dNTP_PopulationKinetics(i,: ) = interp1( Tv_dNTP, Cv_dNTP(:,1) ,dNTP_TimeSpan  ,  'PCHIP');
    DSB_PopulationKinetics(i,: )  = interp1( Tv_DSB ,  Cv_DSB  ,DSB_TimeSpan  ,  'PCHIP');
    ParaPopulation_permutated(i,:) =  Para_permutated;
end

save([ outdir 'Population_', num2str(Casenum),'.mat'], 'ParaPopulation_permutated','TumorVolume_PopulationKinetics','Anabolites_PopulationKinetics','DSB_PopulationKinetics');

end
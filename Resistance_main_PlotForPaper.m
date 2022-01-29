clear
close all
addpath('../DataSet_TGI')

%Parameters fall into four category
%The one element means the category is chosen.
TSinhibion_flag_col = [0,1,1,1,1,0,0,0,0];
Anabolites_flag_col = [1,0,0,0,0,1,0,0,0];
dNTP_flag_col =         [0,0,0,0,0,0,0,1,1];
DSB_flag_col =           [0,0,0,0,1,1,1,0,1];
TSInhibition_TestedPara_col = {[], ["alpha", "59","95","09"],["G_0","cat","08",  "59","95","09"], [ "59","95","09"],...
   ["59","95","09"], [],[],[],[]};
Anabolites_TestedPara_col = {["V_max_54","K_m_54"],[],[],[],[],["V_max_54","K_m_54"],[],[],[]};
ParaPopulation_size =  800; %total number of parameters that are sampled out.
%A = figure('Position', get(0, 'Screensize'));
Case_num = 3;%number of cases per batch
%batch_num = 2;
for i = 1:9
    fprintf('case: %d', i)
    %count = i-3*(batch_num-1);  %for plot,starting from 1 for each batch
    TSInhibition_flag = TSinhibion_flag_col(i);
    Anabolites_flag = Anabolites_flag_col(i);
    dNTP_flag = dNTP_flag_col(i);
    DSB_flag = DSB_flag_col(i); 
    TSInhibition_TestedPara  = TSInhibition_TestedPara_col(i);
    TSInhibition_TestedPara  = [TSInhibition_TestedPara{:}];
    Anabolites_TestedPara = Anabolites_TestedPara_col(i);
    Anabolites_TestedPara  = [Anabolites_TestedPara{:}]; 
    ResistanceAnalysis(Anabolites_flag, Anabolites_TestedPara,...
        TSInhibition_flag ,TSInhibition_TestedPara,...
        dNTP_flag,DSB_flag, ParaPopulation_size )
  
end
%fullpath = strcat(folderpath, "\Allcase.png" );
%saveas(A,fullpath )
newdir = fullfile('..',  'plot' , 'Resistance' ) ;   % Your destination folder
%newdir  =strcat(FolderName , ['/', 'q', num2str(idx)  ]   );
if ~isfolder( newdir ) 
   mkdir( newdir  ) ;
else
   delete(  fullfile(newdir  ,'*' )   )
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  saveas( FigHandle, fullfile(newdir , ['Case_', FigName '.png']));
end


function ResistanceAnalysis(Anabolites_flag, Anabolites_TestedPara,...
    TSInhibition_flag ,TSInhibition_TestedPara,...
    dNTP_flag,DSB_flag, ParaPopulation_size ) 

UnitConversion = [  10^6/130.077, 24*60,10^-3 ,7*24*60   ];
UnitConversion_cell = num2cell(UnitConversion  );
[Plasma_UnitConversion, day2min  ,mm3tocm3  ,week2min  ] = UnitConversion_cell{:};
Table_dir = '../DataSet_all.xlsx';
[T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir) ;
t_dNTP = cell2mat(T_allsheets_output_cell(9) ) ;
t_TS = cell2mat(T_allsheets_output_cell(6) ) ;
t_DSB   = cell2mat(T_allsheets_output_cell(10) ) ; %last data point
t_anabolites = t_TS ;
Data_array_Tv = importdata('50mgkg_EveryTheOtherDay_LS174T_resistant.mat');
Dose_info = [15*60*24  5  48*60  93/2*Plasma_UnitConversion];
theta_Tv_50_sensitive =[   893.53059   50.051548422 0.835170629 0.040987103   0.000269765  5.826570136 0.000019946];
theta_Tv_50_resistance = [  893.53059  69.380333945 0.526513177 49.468042148 0.000249874 7.717727070 0.000004283];

Tv_initial = 125*mm3tocm3;
t_Tumor = Data_array_Tv(:,1); %The prediction doesn't proceed beyond the experimental data

theta_FNUC = [25.17731647	2321.536075];
V_max_54 = theta_FNUC (1);
K_m_54 = theta_FNUC (2);
theta_TS = [21.94138832	;0.198449665;17.39699251;0.032339413; 0.039127057;2.886080266;0.17527698;2.0212;1.0345;0.0186]; 
theta_TS_cell = num2cell(theta_TS);
[K_dUMP, G_0, k_cat,k_08,k_95, k_59 ,k_09, alpha, k_d,TS_0] = theta_TS_cell{:}; 

theta_dNTP = [
       0.21351;32.50090;3.71407;73.34474;3.77027;0.30262;3.33454;1.35311;0.27614;35.66968;6.04530;0.07853;  0.15
        ];
theta_dNTP_cell = num2cell(theta_dNTP);
[k1, k2,k3,k4,k5,k6,k7,k8,k9, k_B,k_A,k10 ,gamma_dNTP ] = theta_dNTP_cell{:};
theta_DSB_delay = 2577.17019; 
lag_DSB   = theta_DSB_delay;
alpha  =2.0212;
Concerned_Para_original= [ V_max_54,K_m_54,G_0, k_cat,k_08,k_95,k_59 ,k_09, k5,lag_DSB,alpha      ] ;
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
dNTP_TimeSpan  = min(t_dNTP):20: max(t_dNTP)+Elongation;
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

ChangedPara = [];
if Anabolites_flag ==1
    ChangedPara = [ChangedPara  Anabolites_TestedPara  ];
    if find(Anabolites_TestedPara == 'V_max_54' )
       LowerBound_fold  = 2;
       ub  = V_max_54; 
       lb = 1/LowerBound_fold*V_max_54;
       Population_Vm54 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
       Population_Vm54= V_max_54*ones(ParaPopulation_size ,1 );
    end
    if find(Anabolites_TestedPara =='K_m_54')
        UpperBound_fold = 2;
        ub = UpperBound_fold* K_m_54;
        lb =  K_m_54;
        Population_Km54 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_Km54 =   K_m_54 *ones(ParaPopulation_size ,1 );
    end
else
     SplitOnCol = num2cell( ones(ParaPopulation_size,1 ) * [ V_max_54,  K_m_54  ],1);%1 :split on column
     [ Population_Vm54,Population_Km54] = SplitOnCol{:}; 
end

if dNTP_flag == 1 % k5 increse the rate of recovery  
    ChangedPara = [ChangedPara "k5"  ];
    UpperBound_fold = 10 ;
    lb =  k5 ;
    ub =UpperBound_fold*k5 ;
    Population_k5 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
else
    Population_k5 =  k5*ones(ParaPopulation_size ,1 ); %vector
end
if TSInhibition_flag ==1 %"G_0","cat","08"  "59","95","09"
    ChangedPara = [ChangedPara TSInhibition_TestedPara ];
    if find(TSInhibition_TestedPara == 'G_0' )
        UpperBound_degree = 10;
        lb =   G_0;
        ub = UpperBound_degree *G_0;
        Population_G0 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_G0 =  G_0*ones(ParaPopulation_size ,1 );
    end 
    if find(TSInhibition_TestedPara =='cat')
        LowerBound_fold = 10 ;
        lb =  1/LowerBound_fold*k_cat;
        ub = LowerBound_fold * k_cat;
        Population_cat =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_cat =  k_cat*ones(ParaPopulation_size ,1 );
    end
    if find(TSInhibition_TestedPara =='08') %dUMP accumulation
        LowerBound_fold = 10;
        ub =   k_08;
        lb = 1/LowerBound_fold *k_08;
        Population_08 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_08 = k_08*ones(ParaPopulation_size ,1 );
    end
    if find(TSInhibition_TestedPara =='59')%dissociation
        UpperBound_fold = 10 ;
        lb =  k_59;
        ub =UpperBound_fold*k_59;
        Population_59 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_59 =  k_59*ones(ParaPopulation_size ,1 );
    end
    if find(TSInhibition_TestedPara =='95') %association
        LowerBound_fold = 10;
        lb = 1/LowerBound_fold*k_95;
        ub = k_95;
        Population_95 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_95 =  k_95*ones(ParaPopulation_size ,1 );%retain value
    end
    if find(TSInhibition_TestedPara =='09') %degradation/elimination 
        UpperBound_fold = 10;
        ub = UpperBound_fold* k_09;
        lb =  k_09;
        Population_09 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_09 =  k_09*ones(ParaPopulation_size ,1 );%retain value
    end
    
    if find(TSInhibition_TestedPara == 'alpha' )
        UpperBound_degree = 10;
        lb =  alpha;
        ub = UpperBound_degree *alpha;
        Population_alpha =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_alpha =  alpha*ones(ParaPopulation_size ,1 );
    end
else
    %if entire TS mechanism is not of interest
   SplitOnCol = num2cell( ones(ParaPopulation_size,1 ) * [ G_0, k_cat,k_08,k_95,k_59 ,k_09,alpha   ] ,1);%1 :split on column
   [ Population_G0, Population_cat,Population_08 ,Population_95 , Population_59, Population_09,  Population_alpha]  = SplitOnCol{:};
end
if DSB_flag == 1 
    ChangedPara = [ChangedPara  "lagDSB"  ];
    UpperBound_fold = 10 ;
    lb =lag_DSB   ;
    ub =UpperBound_fold*lag_DSB   ;
    Population_Tdsb =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
else
    Population_Tdsb =  lag_DSB  *ones(ParaPopulation_size ,1 ); %vector
end

ParaPopulation_unpermutated = [ Population_Vm54   Population_Km54  Population_G0 Population_cat  Population_08   Population_95...
    Population_59 Population_09  Population_k5  Population_Tdsb , Population_alpha];
%V_max_54,K_m_54,G_0, k_cat,k08,k_95,k_59 ,k_09, k5,lag_DSB  
ParaPopulation_permutated = ParaPopulation_unpermutated;


for i = 1: size(  Concerned_Para_original ,2) %shuffle the parameter population matrix
    premuted_column = randperm( ParaPopulation_size  );
    ParaPopulation_permutated(:,i) =  ParaPopulation_unpermutated( premuted_column ,i   );
end
for i = 1:ParaPopulation_size %all the simulation should be have the same time span
    %disp(num2str(i));
    theta_variation = ParaPopulation_permutated(i,:);
    [T_UPToFreeTS,Cv_UPToFreeTS] = Resistance_kinetics_uptoFreeTS(theta_variation,Dose_info);
    [T_dNTP,Cv_dNTP] = Resistance_kinetics_dNTP_plot(theta_variation,Dose_info);
    [T_DSB,Cv_DSB] = Resistance_kinetics_DSB_plot(theta_variation,Dose_info);
    [T_Tumor,Cv_tumor] = Resistance_kinetics_Tv(theta_variation,theta_Tv_50_sensitive, Tv_initial, Dose_info); 
    Anabolites_PopulationKinetics(i,:) =  interp1( T_UPToFreeTS , Cv_UPToFreeTS(:,5) ,Anabolites_TimeSpan  ,  'PCHIP');
    TSInhibtion_PopulationKinetics(i,:) =  interp1( T_UPToFreeTS , Cv_UPToFreeTS(:,end) ,TSInhibition_TimeSpan  ,  'PCHIP');
    dNTP_PopulationKinetics(i,: ) = interp1( T_dNTP, Cv_dNTP(:,1) ,dNTP_TimeSpan  ,  'PCHIP');
    DSB_PopulationKinetics(i,: )  = interp1( T_DSB ,  Cv_DSB  ,DSB_TimeSpan  ,  'PCHIP');
    TumorVolume_PopulationKinetics(i,:) = interp1(T_Tumor , Cv_tumor   ,TumorVolume_TimeSpan,  'PCHIP');
 
end
%  use all(TumorVolume_PopulationKinetics>=0,'all'  ) to see if all the
%  entries are nonegative.
TumorVolume_PopulationKinetics(TumorVolume_PopulationKinetics<0) = 0;

[T_UPtoFreeTS,Cv_UPToFreeTS] = Resistance_kinetics_uptoFreeTS(Concerned_Para_original,Dose_info);
[T_dNTP,Cv_dNTP] = Resistance_kinetics_dNTP_plot(Concerned_Para_original,Dose_info);
[T_DSB,Cv_DSB] = Resistance_kinetics_DSB_plot(Concerned_Para_original,Dose_info);
Anabolites_original  = interp1( T_UPtoFreeTS ,Cv_UPToFreeTS(:,5) , Anabolites_TimeSpan    ,  'PCHIP'); %percentage free TS
TSInhibition_original  = interp1( T_UPtoFreeTS ,Cv_UPToFreeTS(:,end) , TSInhibition_TimeSpan    ,  'PCHIP'); %percentage free TS
dNTP_original =  interp1( T_dNTP ,Cv_dNTP(:,1) ,dNTP_TimeSpan  ,  'PCHIP');
DSB_original =  interp1( T_DSB,Cv_DSB ,DSB_TimeSpan  ,  'PCHIP');

Concerned_Para_original_round = round(Concerned_Para_original,3);
ParaLabel =[ "V_{max,54}","K_{m,54}", "G_{0,dUMP}","k_{cat,dUMP}","k_{08,dUMP}", "k_{95,complex}","k_{59,complex}","k_{09,complex}","k_{5,dNTP}","lag_{DSB} ","alpha"];%string array
%  Turn nummeric array to string array using string() instead of num2str()!!!!! Using num2str() will convert the entire numeric array into one character.

%----------------------------------
%  The following codes are used to summarize the analyze the simulation results of all the parameter samples. 
%
% `idx` points to the parameter set that can give rise to the simulation result closest to the resistant response. 
%----------------------------------

[CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN ] = Resistance_InferentialAnalysis(...
    TumorVolume_PopulationKinetics, ParaPopulation_size, ParaPopulation_permutated    );%CI: confidence interval
[T_tumor_50_resistance,Cv_tumor_50_resistance] = Resistance_kinetics_Tv(Concerned_Para_original,theta_Tv_50_resistance, Tv_initial, Dose_info_Tv);
[T_tumor_50_sensitive,Cv_tumor_50_sensitive] = Resistance_kinetics_Tv(Concerned_Para_original,theta_Tv_50_sensitive , Tv_initial, Dose_info_Tv);
[idx , ClosetTimeCourse] = CompareTimeCourse(T_tumor_50_resistance,Cv_tumor_50_resistance,...
    TumorVolume_PopulationKinetics ,TumorVolume_TimeSpan );
para_ClosetTimeCourse =  ParaPopulation_permutated(idx,:);
[T_tumor_close,Cv_tumor_50_close] = Resistance_kinetics_Tv(para_ClosetTimeCourse,theta_Tv_50_sensitive, Tv_initial, Dose_info_Tv);
 
%fprintf("Case: %d\n", count+3*(batch_num-1))
%subplot(Case_num,3,3*(count-1)+1) 
figure('Position',[1,399,1440,398])
subplot(1,3,1) 
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
grid off ;box off
xlabel('Time(week)')
ylabel('fold change')
xticks( 0:week2min: fix( Dose_info(1) /week2min)*week2min )
xticklabels( (  num2cell( 0:1: fix( Dose_info(1) /week2min) )) )
% title("Tumor Volume kinetics"+newline + "original parameters: "+ OriginPara_Display )  
disp('Tumor volume')
fprintf("\t Max: %s \n",  strjoin( max_legend) )
fprintf("\t Min: %s \n",  strjoin( min_legend) )
fprintf("\t Closest: %s \n",  strjoin( Closest_legend) )
% if  rem(count,3) == 1
%     legend([f1,f2,f3,f4,ha] , {"Simiulated response of resistant LS174T xenograft tumor" ,"Simiulated response of sensitive LS174T xenograft tumor",...
%         'Mean time course',...
%         'Time course closest to the actual resistent response',...
%         ' ', '90% confidence interval'},...
%         'Location','northwestoutside','FontSize', 11)
% end
  
[T_UPtoFreeTS_closeTv,Cv_UPToFreeTS_closeTv] = Resistance_kinetics_uptoFreeTS(para_ClosetTimeCourse,Dose_info);
[T_DSB_closeTv,Cv_DSB_closeTv] = Resistance_kinetics_DSB_plot(para_ClosetTimeCourse,Dose_info);
Anabolites_closeTv  = interp1(T_UPtoFreeTS_closeTv,Cv_UPToFreeTS_closeTv(:,5) , Anabolites_TimeSpan    ,  'PCHIP'); %percentage free TS
DSB_closeTv =  interp1( T_DSB_closeTv,Cv_DSB_closeTv ,DSB_TimeSpan  ,  'PCHIP');
%para_MAX and para_MIN are column vectors
[CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN ] = Resistance_InferentialAnalysis(...
    Anabolites_PopulationKinetics, ParaPopulation_size,  ParaPopulation_permutated     );%CI: confidence interval
[T_data_xtick, T_data_xlabel] = XtickLabelCreation( t_anabolites  , t_anabolites /60 ); %T_data , T_data/24/60

%subplot(Case_num,3,3*(count-1)+2)
subplot(1,3,2) 
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
grid off ;box off
xlabel('Time(h)')
ylabel('pmol/mg')
xticks( T_data_xtick )
xticklabels(T_data_xlabel )
disp('Anabolites')
fprintf("\t Max: %s \n",  strjoin( max_legend) )
fprintf("\t Min: %s \n",  strjoin( min_legend) )
fprintf("\t Closest: %s \n",  strjoin( Closest_legend) )


%subplot(Case_num,3,3*(count-1)+3)
subplot(1,3,3)
[T_data_xtick, T_data_xlabel] = XtickLabelCreation(t_DSB ,t_DSB /60 ); %T_data , T_data/24/60
[CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN ] = Resistance_InferentialAnalysis(...
    DSB_PopulationKinetics, ParaPopulation_size,   ParaPopulation_permutated     );%CI: confidence interval
if C_MAX == C_MIN
    plot(  DSB_TimeSpan, DSB_original  , 'LineWidth',2,'DisplayName' ,'original'   );
else
    f1 = plot( DSB_TimeSpan, DSB_original ,'b', 'LineWidth',2    );
    hold on
    f = plot(  DSB_TimeSpan,C_mean,'--','color',[0.93,0.69,0.13], 'LineWidth',2    ); 
    hold on
    f2 =  plot( DSB_TimeSpan, DSB_closeTv,':' ,'color',[0.49,0.18,0.56], 'LineWidth',2    ); 
    [ha hb hc] = shadedplot( DSB_TimeSpan' ,transpose(CI(:,1)),transpose(CI(:,2)), 'c','none');
    hold off
    max_legend   = CompareWithOriginalValue( para_MAX, Concerned_Para_original,ParaLabel);
    min_legend  = CompareWithOriginalValue( para_MIN, Concerned_Para_original,ParaLabel);

end
grid off ;box off

if DSB_flag ==1
    xlabel('Time(week)')
    ylabel('count(thousands)')
    xticks( 0:week2min: fix( Dose_info(1) /week2min)*week2min )
    xticklabels( (  num2cell( 0:1: fix( Dose_info(1) /week2min) )) )
else
    xlabel('Time(day)')
    ylabel('count(thousands)')
    xticks( min(t_DSB) :24*60: max(t_DSB) )
    xticklabels( num2cell(  min(t_DSB) : 1  : max(t_DSB)/60/24 )) 
end
disp('DSB')
fprintf("\t Max: %s \n",  strjoin( max_legend) )
fprintf("\t Min: %s \n",  strjoin( min_legend) )
fprintf("\t Closest: %s \n",  strjoin( Closest_legend) )
  
end



 

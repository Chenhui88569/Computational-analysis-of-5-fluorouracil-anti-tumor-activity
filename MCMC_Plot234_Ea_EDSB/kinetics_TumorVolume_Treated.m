function Ysim = kinetics_TumorVolume_Treated(theta,Model_Index)
% c_Tv_table= readtable('.\TV_100mgkgweekly.xls'); %relative tumor volume(fold change)
% c_Tv_array = table2array(c_Tv_table);
%num_transit =4 ;
day2min = 24*60;
% mm3tocm3 = 10^-3;
% Anabolites_UnitConversion = 10^6/130.077;
%Dose = 20;
gamma = 0.2;
gamma_hill = gamma;
gamma_DSB = gamma;

i = Model_Index;
AllModel_info   = CollectDataIntoNestedStructure ;
Data_source = AllModel_info(i).DataSource;
Data_array = importdata( Data_source) ;
if Model_Index== 1
    t_data_TV = Data_array(:,1) * day2min ;
else
    t_data_TV =    Data_array(:,1);
end
Dose_info =  AllModel_info(i).Dose_regime;
InitialTV  =AllModel_info(i).InitialTv;
UpstreamKinetics = AllModel_info(i).UpstreamKinetics;
load(UpstreamKinetics, 'T_sim_treated', 'Cv_anabolites_treated','Cv_DSB_treated','T_sim_control', 'Cv_anabolites_control','Cv_DSB_control' )

% if Dose ==100
%     InitialTV  = 125*mm3tocm3;
%     Data_array = importdata('DataSet_TGI/TV_100mgkgweekly.mat');
%     Dose_info =  [45*day2min	3	7*day2min	 92.9936306]; % In line with data
%     t_data_TV = Data_array(:,1) * day2min ; 
%     load Time_Anabolites_DSB_100_weekly_control.mat  T_sim   Cv_anabolites  Cv_DSB
%     T_sim_control = T_sim;
%     Cv_anabolites_control =  Cv_anabolites;
%     Cv_DSB_control  = Cv_DSB;
%     load Time_Anabolites_DSB_100_weekly_treated.mat  T_sim   Cv_anabolites  Cv_DSB
%     
% elseif Dose == 50 % need Data_array, Dose_info
%     i = 2;
%     Data_source = AllModel_info(i).DataSource;
%     Data_array = importdata( Data_source) ;
%     t_data_TV =    Data_array(:,1);
%     Dose_info =  AllModel_info(i).Dose_regime;
%     InitialTV  =AllModel_info(i).InitialTv;
%     load Time_Anabolites_DSB_50_daily_control.mat  T_sim   Cv_anabolites  Cv_DSB
%     T_sim_control = T_sim;
%     Cv_anabolites_control =  Cv_anabolites;
%     Cv_DSB_control  = Cv_DSB;
%     load Time_Anabolites_DSB_50_daily_treated.mat  T_sim   Cv_anabolites  Cv_DSB
% 
%     %Data_array = importdata('50mgkg_daily.mat');
% %     Dose_info = [ Data_array(end,1)	  3     day2min  93/2*Anabolites_UnitConversion]; %weekly
% %     Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_50Daily.mat'); %to save time; take advantage of the fixed duration
%     %Data_array = load('50mgkg_EveryTheOtherDay_LS174T_resistant.mat');
%     %Dose_info = [ 15*60*24  5  48*60  93/2*Anabolites_UnitConversion]; %weekly
% 
% elseif Dose ==20
% %     InitialTV  =100 *mm3tocm3  ;  %100 *mm3tocm3 ,  0.0565
% % %     Data_array = importdata('.\20mgkg_TwiceAweek_Data1_SW620.mat' );%one sheet, one cell
% % %     Dose_info = [    3*7*day2min+5*day2min	3*2-1	7*day2min/2	  93/5*Anabolites_UnitConversion];
% % %     Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_20TwiceAweek.mat');
% %    Data_array = load('20mgkg_ThreeTimesAweek_HT29.mat');%one sheet, one cell 
% %    Dose_info = [   3*7*day2min	 3*3-1	7*day2min/3 	 93/5*Anabolites_UnitConversion]; 
% % %    Data_array = importdata('20mgkg_EveryTheOtherDay_SW620.mat');%one sheet, one cell
% % %    Dose_info = [   2*7*day2min    6     2*day2min  	 93/5*Anabolites_UnitConversion];
% % %    Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_20_EveryTheOtherDay.mat'); %to save time; take advantage of the fixed duration
%     i = 3;
%     Data_source = AllModel_info(i).DataSource;
%     Data_array = importdata( Data_source) ;
%     t_data_TV =    Data_array(:,1);
%     Dose_info =  AllModel_info(i).Dose_regime;
%     InitialTV  =AllModel_info(i).InitialTv;
%     load Time_Anabolites_DSB_20_TwiceAWeek_4Weeks_control.mat   T_sim   Cv_anabolites  Cv_DSB
%     T_sim_control = T_sim;
%     Cv_anabolites_control =  Cv_anabolites;
%     Cv_DSB_control  = Cv_DSB;
%     load Time_Anabolites_DSB_20_TwiceAWeek_4Weeks_treated.mat  T_sim   Cv_anabolites  Cv_DSB
% 
% 
% 
% 
% elseif Dose ==5
%    InitialTV  = 56.76*mm3tocm3;  
%    %InitialTV_NoTreatment  = 46.85*mm3tocm3;  
%    Data_array = load('5mgkg_daily.mat');%one sheet, one cell 
%    Dose_info = [       9*day2min 	8	 day2min 	 92.9936306*Anabolites_UnitConversion/20]; 
%  
% 
% end
Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:};
AnabolitesDSB_TimeSpan = T_sim_treated;
Anabolites_Course = Cv_anabolites_treated;
DSB_Course = Cv_DSB_treated;
AnabolitesDSB_TimeSpan_NoTreatement = T_sim_control;
DSB_Course_NoTreatement  =  Cv_DSB_control;

NonNaN_Rowindex = find( all(~isnan(Data_array), 2) ); 
Relative_data_Treatment = Data_array(:,2);
Relative_data_NoTreatment = Data_array(NonNaN_Rowindex,3);
    
c0_Tv=[ InitialTV      0    InitialTV  ]      ;%initial condition

Cv_Pcell  = [];
T_Pcell  = [];
count = 1;

if ismember(Model_Index, [4,6,7])
    theta_effect  = [893.420616623237;0.00380695282947939;0.172083941222444;0.245506067910982];
    theta_ControlGroupPara = theta(1:3);
else
    theta_effect = theta(1:4);
    theta_ControlGroupPara = theta(5:end);
end
[Tv,Cv]=ode15s(@fun_tumor_apop, t_data_TV , c0_Tv ); %time span
Cv_ProfPlusNonprof   = sum(Cv(:,[1,2]),2);
Cv_ProfPlusNonprof  = real(Cv_ProfPlusNonprof);
Cv_ControlGroup  =   real(Cv(NonNaN_Rowindex,3));

% Ysim =    { Cv_ProfPlusNonprof./c0_Tv(1)  - Relative_data_Treatment  ;
%  Cv_ControlGroup./c0_Tv(3)  - Relative_data_NoTreatment}  ;  %best
Ysim =    { Cv_ProfPlusNonprof -  Relative_data_Treatment* InitialTV
 Cv_ControlGroup - Relative_data_NoTreatment * InitialTV }  ;  %best

    
    function dC = fun_tumor_apop(t_ode,c)
    Cv_Pcell(count) = c(1);
    T_Pcell(count) = t_ode;
    [T_Pcell_unique, ia, ic] = unique(T_Pcell,'sorted');
    Cv_Pcell_2 =  Cv_Pcell(ia);
    
    %%"T_{tumor}","E_{max}","EC_{50}", "\gamma_{DSB}","k_{death}", 

    dcdt = zeros(3,1);
    lag_tumor = theta_effect(1);
    IC50 = theta_effect(2);
    E_max_damage = theta_effect(3);
    EC_50 =  theta_effect(4);
    lambda_g =  theta_ControlGroupPara(1)  ;         %         theta_ControlGroupPara(1)/day2min  ;
    P_max =   theta_ControlGroupPara(2) ;  %theta_ControlGroupPara(2) ; %  1.41650 ,    carrying capacity 
    lambda_d =  theta_ControlGroupPara(3) ; % theta_ControlGroupPara(3)/day2min;%0.0432/60/24 ;  
 
    
    DSB_q_Treatment = interp1(AnabolitesDSB_TimeSpan, DSB_Course , t_ode, 'PCHIP');
    DSB_q_NoTreatment = interp1(AnabolitesDSB_TimeSpan_NoTreatement, DSB_Course_NoTreatement , t_ode, 'PCHIP');

    c_Anabolites_q = interp1( AnabolitesDSB_TimeSpan, Anabolites_Course , t_ode, 'PCHIP');
    DSB_q_Treatment_lag  = DetermineLagingTime(t_ode,DoseFrequency, lag_tumor ,AnabolitesDSB_TimeSpan, DSB_Course);
    DSB_q_NoTreatment_lag  = DetermineLagingTime(t_ode,DoseFrequency, lag_tumor ,AnabolitesDSB_TimeSpan_NoTreatement, DSB_Course_NoTreatement);
    if  DSB_q_NoTreatment_lag == 0
        DSB_q_Treatment_lag = DSB_Course(1);
        DSB_q_NoTreatment_lag = DSB_Course(1);
    end    
    c_Prof_lag_q = DetermineLagingTime(t_ode,DoseFrequency, lag_tumor ,T_Pcell_unique,  Cv_Pcell_2);
     
    Growth_kinetics_Treatment =  lambda_g*c(1)*(1- c(1)/P_max);
    Growth_kinetics_NoTreatment =  lambda_g*c(3)*(1- c(3)/P_max);
   
    E_deviation =@(DSB_0, DSB) ( DSB -DSB_0 )/DSB_0  ;
    E_damage =@(DSB_deviation)   E_max_damage*DSB_deviation^gamma_DSB/(EC_50^gamma_DSB + DSB_deviation^gamma_DSB    );
    %E_damage = @(DSB) c_DSB_TimeCourse(1) /DSB  /( mu_DSB+ c_DSB_TimeCourse(1) /DSB   );
    %E_drug = @(c_Anabolites)  E_max.* c_Anabolites^gamma_hill /(EC50^ gamma_hill+c_Anabolites^gamma_hill  ) ;
    E_drug = @(c_Anabolites)  IC50^gamma_hill /( IC50^gamma_hill +c_Anabolites^gamma_hill) ; %exposure-effect
 
    dcdt(1) = E_drug(c_Anabolites_q) *Growth_kinetics_Treatment -  (1+ E_damage ( E_deviation(DSB_q_NoTreatment, DSB_q_Treatment) ) )  * lambda_d*c(1) ; %proliferating cells 
    dcdt(2) = E_damage ( E_deviation(DSB_q_NoTreatment, DSB_q_Treatment) )* lambda_d*c(1)  -   E_damage (E_deviation(DSB_q_NoTreatment_lag, DSB_q_Treatment_lag))* lambda_d*c_Prof_lag_q;
    dcdt(3) = Growth_kinetics_NoTreatment - lambda_d*c(3) ;
    count = count+1;
    dC=dcdt;
    end
end 
 
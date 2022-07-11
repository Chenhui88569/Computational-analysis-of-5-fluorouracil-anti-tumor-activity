function obj_output  = kinetics_TumorVolume_estimation(theta)
 addpath('../TimeCourses_DSB_Anabolites')
% c_Tv_table= readtable('.\TV_100mgkgweekly.xls'); %relative tumor volume(fold change)
% c_Tv_array = table2array(c_Tv_table);
%num_transit =4 ;
day2min = 24*60;
mm3tocm3 = 10^-3;
Anabolites_UnitConversion = 10^6/130.077;
Dose = 20;
gamma = 0.2;
gamma_hill =gamma;
gamma_DSB =gamma;

%for estimating 20 three times week 
% theta_Tv(1:4) = [893.53059  0.00445 0.24860 0.26015];
% theta_ControlGroupPara  = theta ;
%for model fitting to 5mg/kg daily
% theta_Tv (1:4) = [893.53059 0.00445 0.24860 0.26015]; 
% theta_ControlGroupPara(1:3) = theta;
 
%for model fitting to 20 every the other day
theta_Tv  = [  893.53059  0.00445 0.24860 0.26015 ];
theta_ControlGroupPara = theta ;

%for model fitting to 50 every the other day
% theta_Tv  = [893.53059    theta];
% theta_ControlGroupPara = [ 0.000249874 7.717727070 0.000004283];
% 
if Dose ==100
    InitialTV  = 125*mm3tocm3;
%     Data_array = importdata('TV_100mgkgweekly.mat');
%     Data_array(:,1)  = Data_array(:,1)*day2min;
    Data_array = importdata('100mgkg_resistance_Colon26.mat');
    Dose_info =  [45*day2min	3	7*day2min	 92.9936306*Anabolites_UnitConversion]; % In line with data
    Time_Anabolites_DSB_array = importdata('Time_Anabolites_DSB_100weekly.mat'); %to save time ;take advantage of the fixed duration
elseif Dose == 50
    InitialTV  = 125*mm3tocm3; 
%     Data_array = importdata('50mgkg_daily.mat');
%     Dose_info = [ Data_array(end,1)	  3     day2min  93/2*Anabolites_UnitConversion]; %weekly
%     Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_50Daily.mat'); %to save time; take advantage of the fixed duration
    Data_array = importdata('50mgkg_EveryTheOtherDay_LS174T_resistant.mat');
    Dose_info = [ 15*60*24  5  48*60  93/2*Anabolites_UnitConversion]; %weekly
    Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_50_EveryTheOtherDay.mat'); %to save time; take advantage of the fixed duration

elseif Dose ==20
    InitialTV  =100 *mm3tocm3  ;  %100 *mm3tocm3 ,  0.0565
%     Data_array = importdata('.\20mgkg_TwiceAweek_Data1_SW620.mat' );%one sheet, one cell
%     Dose_info = [    3*7*day2min+5*day2min	3*2-1	7*day2min/2	  93/5*Anabolites_UnitConversion];
%     Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_20TwiceAweek.mat');
   Data_array = importdata('20mgkg_ThreeTimesAweek_HT29.mat');%one sheet, one cell 
   Dose_info = [   3*7*day2min	 3*3-1	7*day2min/3 	 93/5*Anabolites_UnitConversion]; 
   Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_20ThreeTimesAWeek.mat'); %to save time; take advantage of the fixed duration
%    Data_array = importdata('20mgkg_EveryTheOtherDay_SW620.mat');%one sheet, one cell
%    Dose_info = [   2*7*day2min    6     2*day2min  	 93/5*Anabolites_UnitConversion];
%    Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_20_EveryTheOtherDay.mat'); %to save time; take advantage of the fixed duration
% 

elseif Dose ==5
   InitialTV  = 56.76*mm3tocm3;  
   %InitialTV_NoTreatment  = 46.85*mm3tocm3;  
   Data_array = importdata('5mgkg_daily.mat');%one sheet, one cell 
   Dose_info = [       9*day2min 	8	 day2min 	 92.9936306*Anabolites_UnitConversion/20]; 
   Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_5_daily.mat'); %to save time; take advantage of the fixed duration
 

end
Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:};

AnabolitesDSB_TimeSpan =Time_Anabolites_DSB_array(:,1);
Anabolites_Course = Time_Anabolites_DSB_array(:,2);
DSB_Course = Time_Anabolites_DSB_array(:,3);
Time_Anabolites_DSB_NoTreatment = importdata('Time_Anabolites_DSB_NoTreatment.mat');
AnabolitesDSB_TimeSpan_NoTreatement =Time_Anabolites_DSB_NoTreatment(:,1);
DSB_Course_NoTreatement  =Time_Anabolites_DSB_NoTreatment(:,3);
t_data_TV = Data_array(:,1) ; 

NonNaN_Rowindex = find( all(~isnan(Data_array), 2) ); 
Relative_data_Treatment = Data_array(:,2);
Relative_data_NoTreatment = Data_array(NonNaN_Rowindex,3);
    
c0_Tv=[ InitialTV      0  InitialTV  ]      ;%initial condition

Cv_Pcell  = [];
T_Pcell  = [];
count = 1;
[Tv,Cv]=ode15s(@fun_tumor_apop, t_data_TV , c0_Tv ); %time span
Cv_ProfPlusNonprof   = sum(Cv(:,[1,2]),2);
Cv_ProfPlusNonprof  = real(Cv_ProfPlusNonprof);
Cv_ControlGroup  =   real(Cv(NonNaN_Rowindex,3));

Ysim =    [ Cv_ProfPlusNonprof./c0_Tv(1)  - Relative_data_Treatment  ;
 [Cv_ControlGroup./c0_Tv(3)  - Relative_data_NoTreatment]*1/2 ]  ;  %best

%obj_output = Ysim;
SSR =  sum( Ysim.^2);
%%obj_output = sqrt(SSR);
ObjFunction = @(num_datapoints, SSR, num_para) num_datapoints*log(SSR/num_datapoints)+2*(num_para+1)*num_datapoints...
    /(num_datapoints-num_para-2);
num_datapoints = size(Ysim,1);
num_para = size( theta,2);
obj_output =   ObjFunction(num_datapoints,SSR,num_para);
    
    function dC = fun_tumor_apop(t_ode,c)
    Cv_Pcell(count) = c(1);
    T_Pcell(count) = t_ode;
    [T_Pcell_unique, ia, ic] = unique(T_Pcell,'sorted');
    Cv_Pcell_2 =  Cv_Pcell(ia);
    
    %%"T_{tumor}","E_{max}","EC_{50}", "\gamma_{DSB}","k_{death}", 

    dcdt = zeros(2,1);
    lag_tumor = theta_Tv(1);
    IC50 = theta_Tv(2);
    E_max_damage = theta_Tv(3);
    EC_50 =  theta_Tv(4);
    lambda_g =  theta_ControlGroupPara(1)  ;         %         theta_ControlGroupPara(1)/day2min  ;
    chi_max =   theta_ControlGroupPara(2) ;  %theta_ControlGroupPara(2) ; %  1.41650 ,    carrying capacity 
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
     
    Growth_kinetics_Treatment =  lambda_g*c(1)*(1- c(1)/chi_max);
    Growth_kinetics_NoTreatment =  lambda_g*c(3)*(1- c(3)/chi_max);
   
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
 
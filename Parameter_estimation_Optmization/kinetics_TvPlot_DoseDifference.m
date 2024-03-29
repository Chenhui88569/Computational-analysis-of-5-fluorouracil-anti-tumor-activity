 function [Tsim,Ysim] = kinetics_TvPlot_DoseDifference(theta,TV_initial, Dose_info) 
 addpath('TimeCourses_DSB_Anabolites')
%Usages:
%1. plot 9 models with paramter estimates
%2. simulation of doses and dose regimes

%Duration is a number 

Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:}; 
Time_Anabolites_DSB_array  = kinetics_uptoDSB_alternative(Dose_info);
AnabolitesDSB_TimeSpan =Time_Anabolites_DSB_array(:,1);
Anabolites_TimeCourse = Time_Anabolites_DSB_array(:,2);
DSB_TimeCourse = Time_Anabolites_DSB_array(:,3);

Time_Anabolites_DSB_NoTreatment = importdata('Time_Anabolites_DSB_NoTreatment.mat');
AnabolitesDSB_TimeSpan_NoTreatement =Time_Anabolites_DSB_NoTreatment(:,1);
DSB_TimeCourse_NoTreatement  =Time_Anabolites_DSB_NoTreatment(:,3);

c0_Tv=[ TV_initial   0] ;  
Cv_Pcell  = [];
T_Pcell  = [];
count = 1;
theta_Tv = theta;
gamma = 0.2;
gamma_hill =gamma;
gamma_DSB =gamma;

[Tv,Cv]=ode15s(@fun_tumor_apop,  [0 Duration], c0_Tv ); %time span
Cv_ProfPlusNonprof   = sum(Cv,2);
Cv_ProfPlusNonprof  = real(Cv_ProfPlusNonprof);
Ysim = Cv_ProfPlusNonprof;
Ysim(:,end+1) = real(Cv(:,2));
Tsim =  Tv;
 
function dC = fun_tumor_apop(t_ode,c)
    Cv_Pcell(count) = c(1);
    T_Pcell(count) = t_ode;
    [T_Pcell_unique, ia, ic] = unique(T_Pcell,'sorted');
    Cv_Pcell_2 =  Cv_Pcell(ia);
    %%"T_{tumor}","E_{max}","EC_{50}", "\gamma_{DSB}","k_{out}", 

    dcdt = zeros(2,1);
    lag_tumor = theta_Tv(1);
    IC50 = theta_Tv(2);
    E_max_damage = theta_Tv(3);
    EC_50 = theta_Tv(4);
    lambda_g =  theta_Tv(5)   ;
    chi_max =  theta_Tv(6);  % ; %  1.41650 ,    carrying capacity 
    lambda_d = theta_Tv(7) ;  
     
    DSB_q_Treatment = interp1(AnabolitesDSB_TimeSpan, DSB_TimeCourse , t_ode, 'PCHIP');
    DSB_q_NoTreatment = interp1(AnabolitesDSB_TimeSpan_NoTreatement, DSB_TimeCourse_NoTreatement , t_ode, 'PCHIP');

    c_Anabolites_q = interp1( AnabolitesDSB_TimeSpan, Anabolites_TimeCourse , t_ode, 'PCHIP');
    DSB_q_Treatment_lag  = DetermineLagingTime(t_ode,DoseFrequency, lag_tumor ,AnabolitesDSB_TimeSpan, DSB_TimeCourse);
    DSB_q_NoTreatment_lag  = DetermineLagingTime(t_ode,DoseFrequency, lag_tumor ,AnabolitesDSB_TimeSpan_NoTreatement, DSB_TimeCourse_NoTreatement);
    if  DSB_q_NoTreatment_lag == 0
        DSB_q_Treatment_lag = DSB_TimeCourse(1);
        DSB_q_NoTreatment_lag = DSB_TimeCourse(1);
    end       
    c_Prof_lag_q = DetermineLagingTime(t_ode,DoseFrequency, lag_tumor ,T_Pcell_unique,  Cv_Pcell_2);
 
    Growth_kinetics =  lambda_g*c(1)*(1- c(1)/chi_max);
    E_deviation =@(DSB_0, DSB) ( DSB -DSB_0 )/DSB_0  ;
    E_damage =@(DSB_deviation)   E_max_damage*DSB_deviation^gamma_DSB/(EC_50^gamma_DSB + DSB_deviation^gamma_DSB    );
    E_drug = @(c_Anabolites)  IC50^gamma_hill /( IC50^gamma_hill +c_Anabolites^gamma_hill) ; %exposure-effect

    dcdt(1) = E_drug(c_Anabolites_q) * Growth_kinetics -(1+  E_damage (E_deviation(DSB_q_NoTreatment, DSB_q_Treatment))) *lambda_d *c(1)   ; %proliferating cells 
    dcdt(2) = E_damage ( E_deviation(DSB_q_NoTreatment, DSB_q_Treatment) )*lambda_d*c(1)  -  E_damage (E_deviation(DSB_q_NoTreatment_lag, DSB_q_Treatment_lag))*lambda_d*c_Prof_lag_q;
    count = count+1;
    dC=dcdt;
    end

end

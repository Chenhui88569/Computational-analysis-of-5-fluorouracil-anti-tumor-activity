 function [Tsim,Ysim] = Plot_kinetics_TumorVolume(theta, Model_Index) 
%Usages:
%1. plot 9 models with paramter estimates
%2. simulation of doses and dose regimes

%Duration is a number 
AllModel_info   = CollectDataIntoNestedStructure ;
Dose_info =  AllModel_info( Model_Index).Dose_regime;
TV_initial  =AllModel_info( Model_Index).InitialTv;
UpstreamKinetics = AllModel_info( Model_Index).UpstreamKinetics;
load(UpstreamKinetics, 'T_sim_treated', 'Cv_anabolites_treated','Cv_DSB_treated','T_sim_control', 'Cv_anabolites_control','Cv_DSB_control' )

Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:}; 
% if Dose == 50
%     load Time_Anabolites_DSB_50_daily_treated.mat  T_sim   Cv_anabolites  Cv_DSB
% else
%     [T_sim, Cv_anabolites,Cv_DSB ]   = RepeatedDoses_kinetics_uptoDSBandFreeTS(Model_Index,0); %The input argument is (model_index,Control_flag)
% end
AnabolitesDSB_TimeSpan = T_sim_treated;
Anabolites_TimeCourse = Cv_anabolites_treated;
DSB_TimeCourse = Cv_DSB_treated;

% Dose_info(end) = 0;
% if Dose ==50
%     load Time_Anabolites_DSB_50_daily_control.mat  T_sim   Cv_anabolites  Cv_DSB
% else 
%     [T_sim, Cv_anabolites,Cv_DSB ] = RepeatedDoses_kinetics_uptoDSBandFreeTS(Dose_info);
% end
AnabolitesDSB_TimeSpan_NoTreatement = T_sim_control;
DSB_TimeCourse_NoTreatement  = Cv_DSB_control; 

c0_Tv=[ TV_initial   0 TV_initial] ;  
Cv_Pcell  = [];
T_Pcell  = [];
count = 1;
gamma = 0.2;
gamma_hill =gamma;
gamma_DSB =gamma;
if ismember(Model_Index, [4,6,7])
    %theta_effect  = [893.420616623237;0.00380695282947939;0.172083941222444;0.245506067910982];
    theta_effect = [893.445064075148;0.00376576201660318;0.175114274349874;0.238015233076624];
    theta_ControlGroupPara = theta(1:3);
else
    theta_effect = theta(1:4);
    theta_ControlGroupPara = theta(5:end);
end
[Tsim,Cv]=ode15s(@fun_tumor_apop, 0:50:Duration , c0_Tv ); %time span
Cv_ProfPlusNonprof   = sum(Cv(:,[1,2]),2);
Cv_ProfPlusNonprof  = real(Cv_ProfPlusNonprof);
Cv_ControlGroup  =   real(Cv(:,3));
Ysim = [Cv_ProfPlusNonprof  Cv_ControlGroup ]./TV_initial;

function dC = fun_tumor_apop(t_ode,c)
    dcdt = zeros(3,1);
    Cv_Pcell(count) = c(1);
    T_Pcell(count) = t_ode;
    [T_Pcell_unique, ia, ic] = unique(T_Pcell,'sorted');
    Cv_Pcell_2 =  Cv_Pcell(ia);
    %%"T_{tumor}","E_{max}","EC_{50}", "\gamma_{DSB}","k_{out}", 
  
    lag_tumor = theta_effect(1);
    IC50 = theta_effect(2);
    E_max_damage = theta_effect(3);
    EC_50 =  theta_effect(4);
    lambda_g =  theta_ControlGroupPara(1)  ;         %         theta_ControlGroupPara(1)/day2min  ;
    P_max =   theta_ControlGroupPara(2) ;  %theta_ControlGroupPara(2) ; %  1.41650 ,    carrying capacity 
    lambda_d =  theta_ControlGroupPara(3) ; % theta_ControlGroupPara(3)/day2min;%0.0432/60/24 ;  

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
 
    Growth_kinetics =  lambda_g*c(1)*(1- c(1)/P_max);
    E_deviation =@(DSB_0, DSB) ( DSB -DSB_0 )/DSB_0  ;
    E_damage =@(DSB_deviation)   E_max_damage*DSB_deviation^gamma_DSB/(EC_50^gamma_DSB + DSB_deviation^gamma_DSB    );
    E_drug = @(c_Anabolites)  IC50^gamma_hill /( IC50^gamma_hill +c_Anabolites^gamma_hill) ; %exposure-effect

    dcdt(1) = E_drug(c_Anabolites_q) * Growth_kinetics -(1+  E_damage (E_deviation(DSB_q_NoTreatment, DSB_q_Treatment))) *lambda_d *c(1)   ; %proliferating cells 
    dcdt(2) = E_damage ( E_deviation(DSB_q_NoTreatment, DSB_q_Treatment) )*lambda_d*c(1)  -  E_damage (E_deviation(DSB_q_NoTreatment_lag, DSB_q_Treatment_lag))*lambda_d*c_Prof_lag_q;
    dcdt(3)=  lambda_g*c(3)*(1- c(3)/P_max) - lambda_d*c(3); 
    count = count+1;
    dC=dcdt;
    end

end

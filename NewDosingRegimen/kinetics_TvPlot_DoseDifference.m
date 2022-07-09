function [T_TimeCourse_final,Ysim] = kinetics_TvPlot_DoseDifference(theta_col, Dose_info, TV_initial ,scatter_injection ,Dose_order) 
% theta is cell
dbstop if error
num_dose_PerCycle = size(Dose_info,1);


Control_flag= 0;
[T_sim_treated, Cv_anabolites_treated, Cv_DSB_treated ]  = RepeatedDoses_kinetics_uptoDSBandFreeTS(Dose_info, Control_flag, scatter_injection,Dose_order);

Control_flag= 1;
[T_sim_control, Cv_anabolites_control, Cv_DSB_control]  = RepeatedDoses_kinetics_uptoDSBandFreeTS(Dose_info, Control_flag, scatter_injection,Dose_order);

gamma = 0.2;
gamma_hill =gamma;
gamma_DSB =gamma;
num_injection = length(scatter_injection);
c0 =[ TV_initial   0  TV_initial] ;

Cv_complex_Timecourse  = [];
Cv_complex_Timecourse(1,:) = c0' ;
Tv_TimeCourse =[];
Tv_TimeCourse(1) = 0;

for indi_injection = 1: num_injection
     if indi_injection < num_injection
         Cv_Pcell  = [];
         T_Pcell  = [];
         count = 1;
         indi_dur = [scatter_injection(indi_injection)   scatter_injection(indi_injection+1)] - scatter_injection(indi_injection) ;% Duration starts from 0
         Dose_idx = Dose_order( indi_injection );
         theta = theta_col(Dose_idx ,:);
         theta_effect = theta(1:4);
         theta_ControlGroupPara = theta(5:end);
         
         T_idx_1  = find(  T_sim_treated ==  indi_dur(1)  + scatter_injection(indi_injection)  );
         T_idx_2 = find(  T_sim_treated ==  indi_dur(2)  +  scatter_injection(indi_injection)  );
         AnabolitesDSB_TimeSpan  =   T_sim_treated(T_idx_1:T_idx_2);
         Anabolites_TimeCourse = Cv_anabolites_treated(T_idx_1:T_idx_2);
         DSB_TimeCourse = Cv_DSB_treated(T_idx_1:T_idx_2);
         
         T_idx_1_control  = find( T_sim_control ==  indi_dur(1)  + scatter_injection(indi_injection)  );
         T_idx_2_control = find(  T_sim_control ==  indi_dur(2)  +  scatter_injection(indi_injection)  );

         AnabolitesDSB_TimeSpan_NoTreatement = T_sim_control(   T_idx_1_control :T_idx_2_control);
         DSB_TimeCourse_NoTreatement  = Cv_DSB_control(   T_idx_1_control :T_idx_2_control);
         
         last_time = Tv_TimeCourse(end);
         
         
         [Tv,Cv]=ode15s(@DifEq, indi_dur,c0);
         Tv= Tv';
         Cv(Cv<=0) =0;
         c0 = Cv(end,:); %reinitialize
         Tv_TimeCourse( end+1 : end+size(Tv,2)  )  = Tv + Tv_TimeCourse(end);
         Cv_complex_Timecourse( end+1: end+size(Tv,2) ,:  ) =  Cv;

     elseif    indi_injection == num_injection
         Cv_Pcell  = [];
         T_Pcell  = [];
         count = 1;
         t_interval_final = [0    sum(Dose_info(:,1))  - scatter_injection(end)  ];
         Dose_idx = Dose_order( indi_injection );
         theta = theta_col(Dose_idx ,:);
          theta_effect = theta(1:4);
         theta_ControlGroupPara = theta(5:end);
         AnabolitesDSB_TimeSpan = T_sim_treated;
         T_idx_1  = find(  T_sim_treated ==  scatter_injection(end) );
         T_idx_2 = find(  T_sim_treated ==   sum(Dose_info(:,1))  );   
         AnabolitesDSB_TimeSpan =  T_sim_treated(T_idx_1:T_idx_2);
         Anabolites_TimeCourse = Cv_anabolites_treated(T_idx_1:T_idx_2);
         DSB_TimeCourse = Cv_DSB_treated(T_idx_1:T_idx_2);
         
         T_idx_1_control  = find( T_sim_control == scatter_injection(end) );
         T_idx_2_control = find(  T_sim_control ==  sum(Dose_info(:,1))  );

         AnabolitesDSB_TimeSpan_NoTreatement = T_sim_control(T_idx_1_control:T_idx_2_control);
         DSB_TimeCourse_NoTreatement  = Cv_DSB_control(T_idx_1_control:T_idx_2_control);
         
         last_time = Tv_TimeCourse(end);
         
         [Tv,Cv]=ode15s(@DifEq,t_interval_final,c0); %15s
        Tv = Tv';
        Cv(Cv<=0) =0;
        Tv_TimeCourse( end+1 : end+size(Tv,2)  )  = Tv + Tv_TimeCourse(end); %lump the time courses of interdose intervals
        Cv_complex_Timecourse( end+1: end+size(Tv,2),: ) =  Cv ;
     end
end
[T_TimeCourse_final, ia_final, ~] = unique(Tv_TimeCourse,'sorted'); %horizontal
Cv_complex_Timecourse = Cv_complex_Timecourse( ia_final,:);

Cv_ProfPlusNonprof   = sum(Cv_complex_Timecourse(:,[1,2]),2);
Cv_ProfPlusNonprof  = real(Cv_ProfPlusNonprof);
Cv_ControlGroup  =   real(Cv_complex_Timecourse(:,3));
Ysim = [Cv_ProfPlusNonprof  Cv_ControlGroup ]./TV_initial;

    function dC = DifEq(t_ode,c)        
  
    DoseFrequency = Dose_info(1,3);
    
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

    DSB_q_Treatment = interp1(AnabolitesDSB_TimeSpan, DSB_TimeCourse , t_ode+  last_time, 'PCHIP');
    DSB_q_NoTreatment = interp1(AnabolitesDSB_TimeSpan_NoTreatement, DSB_TimeCourse_NoTreatement , t_ode+  last_time, 'PCHIP');

    c_Anabolites_q = interp1( AnabolitesDSB_TimeSpan, Anabolites_TimeCourse , t_ode+  last_time, 'PCHIP');
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


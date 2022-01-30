function [c, ceq] = AUCe_NonlinearConstraint(theta)
%'nonlcon': Function handle corresponding to the nonlinear constraint function. The constraint function must accept a vector x and return two vectors: c, the nonlinear inequality constraints, and ceq, the nonlinear equality constraints. If one of these constraint functions is empty, nonlcon must return [] for that function.
% function [c,ceq] = unitdisk(x)
% c = x(1)^2 + x(2)^2 - 1;
% ceq = []
day2min  = 60*24; 
FinalTime = 7*day2min; %over a week
inter_timespan= 1:50:FinalTime;
Control_file = 'Time_Anabolites_DSB_NoTreatment.mat';
Time_Anabolites_DSB_control = importdata(Control_file); %to save time ;take advantage of the fixed duration
AnabolitesDSB_TimeSpan_control = Time_Anabolites_DSB_control(:,1);
DSB_Course_control = Time_Anabolites_DSB_control(:,3);
gamma_DSB   = 0.2;
gamma_hill = 0.2;
E_deviation =@(DSB_0, DSB) ( DSB -DSB_0 )./DSB_0  ; 

%AUCe_DSB_against_ub = 1.5989e-05; %AUCe of model over a week,  
%AUCe_anabolites_against_ub = 3.4904; 
AUCe_DSB_against_lb = 3.5154e-05; %AUCe of model over a week, compare to the model 3
%AUCe_anabolites_against_lb =2.0441; 
para_individual =[450.009370742   0.00445 0.24860 0.26015     0.000146106   1.666387293   0.000010027]; %model 3
para_individual_cell  = num2cell( para_individual );
Time_Anabolites_DSB_array =  importdata('Time_Anabolites_DSB_20ThreeTimesAWeek.mat'); %to save time; take advantage of the fixed duration
c =zeros(1,1);
ControlGroupRelatedEstimation = 1; %test
if ControlGroupRelatedEstimation==  1
    AUCe_DSB_against_lb = 3.5154e-05; %AUCe of model over a week, compart to the model 3 
    [~, IC50,E_max_damage, EC50,~, ~,~ ]=   para_individual_cell{:};
    E_damage =@(DSB_deviation , lambda_d)   lambda_d.*E_max_damage.*DSB_deviation.^gamma_DSB./(EC50.^gamma_DSB + DSB_deviation.^gamma_DSB    );  
    E_anabolites = @(c_Anabolites)  IC50.^gamma_hill ./( IC50.^gamma_hill +c_Anabolites.^gamma_hill) ; %exposure-effect 
elseif  ControlGroupRelatedEstimation ==  0
    [~, ~,~, ~,~, ~, lambda_d ]=   para_individual_cell{:};
    E_damage =@(DSB_deviation , E_max_damage,EC50)   lambda_d.*E_max_damage.*DSB_deviation.^gamma_DSB./(EC50.^gamma_DSB + DSB_deviation.^gamma_DSB    );
    E_anabolites = @(c_Anabolites,IC50)  IC50.^gamma_hill ./( IC50.^gamma_hill +c_Anabolites.^gamma_hill) ; %exposure-effect 
 end
AnabolitesDSB_TimeSpan =Time_Anabolites_DSB_array(:,1);
Anabolites_Course = Time_Anabolites_DSB_array(:,2);
DSB_Course = Time_Anabolites_DSB_array(:,3);
DSB_Course_q =  interp1( AnabolitesDSB_TimeSpan  ,DSB_Course ,    inter_timespan  ,'pchip'  );
Anabolites_Course_q =  interp1( AnabolitesDSB_TimeSpan  ,Anabolites_Course ,    inter_timespan  ,'pchip'  );
DSB_Course_control_q =  interp1(  AnabolitesDSB_TimeSpan_control  ,DSB_Course_control ,    inter_timespan  ,'pchip'  );
DSB_deviation = E_deviation(DSB_Course_control_q,   DSB_Course_q );
%apply variables in the intermediate functions E_damage and E_anabolites before constructing the nonlinear
%inequality constraints.
if ControlGroupRelatedEstimation==  1
    E_damage_individual = E_damage(DSB_deviation,theta(3) );  
    AUC_E_DSB = trapz(  inter_timespan./day2min , E_damage_individual  );
    c(1)  =  -real(AUC_E_DSB) +  AUCe_DSB_against_lb ;%<0
   % c(2)  =  real(AUC_E_DSB) - AUCe_DSB_against_ub ;%<0
elseif ControlGroupRelatedEstimation ==  0 %estimate dose related para
    E_damage_individual = E_damage(DSB_deviation,theta(2),theta(3) );
    E_Anabolites_individual =   E_anabolites( Anabolites_Course_q, theta(1)); 
    AUC_E_DSB = trapz(  inter_timespan./day2min , E_damage_individual  );
    AUC_E_A_complement = trapz(  inter_timespan./day2min , E_Anabolites_individual  );
    AUC_E_A =  1*7 -  AUC_E_A_complement ;
    c(1) =  -real(AUC_E_DSB) +AUCe_DSB_against_ub ;
    c(2) =  -real(AUC_E_A)+AUCe_anabolites_against  ;
end
 
ceq = [];

end


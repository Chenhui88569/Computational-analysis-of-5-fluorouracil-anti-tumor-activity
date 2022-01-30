function COL  = kinetics_uptoDSB_alternative(Dose_info)
% Dose_info = [
%     45*day2min	3	7*day2min	 93/20*plasma_UnitConversion
%     45*day2min	3	7*day2min	 93/10*plasma_UnitConversion
%     45*day2min	3	7*day2min	 93/5*plasma_UnitConversion
%     45*day2min	3	7*day2min	 93/2*plasma_UnitConversion
%     45*day2min	3	7*day2min	 93*plasma_UnitConversion
%     45*day2min	3	7*day2min	 0
%     ];%Duration, Interval_TotalNum,DoseFrequency ,Dose
% %Dose_info (:,1:3)  = ones(6,1)* [1*day2min  0   7*day2min	 ];
Dose_info_focus = Dose_info;
 [T_UPtoFreeTS,Cv_UPtoFreeTS]=kinetics_uptoFreeTS(Dose_info_focus); %time span
 theta_DSB = [ 2577.17019
10.50768
151.94553
10.06779
194.02545
0.12550
2.79398];
[T_DSB,Cv_DSB] = kinetics_DSB_plot(theta_DSB,Dose_info_focus); 

Cv_DSB_q= interp1(T_DSB , Cv_DSB , T_UPtoFreeTS,  'PCHIP');
COL = [ T_UPtoFreeTS,    Cv_UPtoFreeTS(:,5)   Cv_DSB_q]; %make sure that (end-1)th column is DSB
%save Time_Anabolites_DSB_20_EveryTheOtherDay.mat  COL 
% X= importdata('Time_Anabolites_DSB_20_EveryTheOtherDay.mat');
% plot(X(:,1),X(:,2));
end

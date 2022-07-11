function  c_q  = DetermineLagingTime(...
    t_ode, DoseFrequency, lag ,T_unique,Cv_2)
%DETERMINELAGINGTIME  generate the circle, periofic curve
if t_ode <= DoseFrequency*1 % During the first interval, the history is zero
    if  t_ode <= lag
        c_q = 0;
    elseif t_ode >= lag
        c_q = interp1(T_unique  ,  Cv_2  , t_ode- lag,  'PCHIP');
    end
elseif   DoseFrequency*1 < t_ode  %subsequent dose
        c_q = interp1(T_unique  ,  Cv_2  , t_ode- lag,  'PCHIP');
        
% elseif   DoseFrequency*1 < t_ode <= DoseFrequency*( Interval_TotalNum-1) %intermediate dose
%     lag_AfterFirstDose = lag + DoseFrequency*fix(t_ode/DoseFrequency) ; %generate the circle, periofic curve
%     c_q = interp1(T_unique  ,  Cv_2  , t_ode- lag_AfterFirstDose,  'PCHIP');
% elseif   DoseFrequency*( Interval_TotalNum-1) < t_ode %last dose
%     lag_AfterFirstDose = lag + DoseFrequency*( Interval_TotalNum-1) ; 
%     c_q = interp1(T_unique  ,  Cv_2  , t_ode- lag_AfterFirstDose,  'PCHIP');
end

end


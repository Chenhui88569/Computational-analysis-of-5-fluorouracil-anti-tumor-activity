function  [Tv_output, Cv_output]  = EstimateTumorGrowth_plot(theta,T_time ,c0_Tv)
[Tv,Cv]=ode15s(@fun_TV,T_time , c0_Tv ); %time span
Tv_output = Tv;
Cv_output =Cv;

    function diff = fun_TV(t,c)
      %  DSB_q = interp1(PlasmaDSB_TimePoint, DSB_TimeCourse , t, 'PCHIP');
        dcdt = zeros(1,1);
        lambda_g = theta(1); 
        P_max = theta(2);
        lambda_d = theta(3);
        dcdt(1)= lambda_g*c(1)*(1- c(1)/P_max) - lambda_d*c(1);        
        diff = dcdt;
    end

end


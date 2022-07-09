function [T_DSB,Cv_DSB,T_uptodNTP,Cv_dNTP] = Resistance_kinetics_DSB_plot(theta_cellular, theta_dNTP, theta_DSB, theta_variation,Dose_info)

theta = theta_DSB;
c0=[
    42.755 % DSB
    ];%initial condition
theta_resistance =   theta_variation;
Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:}; 
%theta_dNTP = [0.21351;32.50090;3.71407;73.34474;3.77027;0.30262;3.33454;1.35311;0.27614;35.66968;6.04530;0.07853];
[T_uptodNTP,Cv_dNTP]= Resistance_kinetics_dNTP_plot(theta_cellular, theta_dNTP, theta_variation,Dose_info);
dNTP_TimeCourse = Cv_dNTP(:,1);
[T_DSB,Cv_DSB]=ode15s(@fun_DSB,[0 Duration], c0 ); %time span
Cv_DSB = real(Cv_DSB);


    function dC = fun_DSB(time_ode,c)
        dcdt = zeros(1,1);
       % theta = [2577.17019 ;10.50768; 151.94553; 10.06779; 194.02545 ;0.12550 ;2.79398 ];
        lag_DSB  = theta(1) ; 
        V_dNTP = theta(2);
        K_dNTP =  theta(3);
        V_HRR =  theta(4);% theta(3);
        K_HRR = theta(5) ;%theta(4);
        k_i =  theta(6);
        k_0 =   theta(7); %theta(6);
        gamma_1  = 0.6;
        theta_resistance_cell = num2cell(theta_resistance);
       [ V_max_54,K_m_54,G_0, k_cat,k_08,k_95,k_59 ,k_09, k5,lag_DSB,alpha, V_HRR, K_HRR,  lag_tumor,E_max_damage, EC_50 ]  = theta_resistance_cell{:};
        lag_DSB   = lag_DSB  *24*60;
       if time_ode <= lag_DSB
            dcdt(1)  = k_0 - ...
                V_HRR.*c(1) /(  K_HRR+c(1) );
        elseif time_ode > lag_DSB
            clag1 = interp1(T_uptodNTP , dNTP_TimeCourse  , time_ode-lag_DSB, 'PCHIP');
             dcdt(1)  = k_0.* (1+ V_dNTP.* clag1^gamma_1 / (  K_dNTP^gamma_1 + clag1^gamma_1  ) )- ...
                V_HRR.*c(1)   /(  K_HRR+c(1) )/(k_i^gamma_1* clag1^gamma_1  + 1  ) ;
        end
        dC=dcdt;
    end

end 

function [T,Cv_DSB,T_UptoFreeTS , C_sim_UpToFreeTS] = RepeatedDoses_kinetics_dNTP_DSB(theta,theta_3to10, Dose_info, Control_flag,scatter_injection,Dose_order)


Cv_dNTP  = [];
T_dNTP = [];


c0=[
    0
    0
    42.755 % DSB
    ];%initial condition

[T_UptoFreeTS , C_sim_UpToFreeTS ] = RepeatedDoses_kinetics_Cellular3to10(theta_3to10, Dose_info, Control_flag, scatter_injection,Dose_order);
PercentagefreeTS = C_sim_UpToFreeTS(:,end);

TotalDuration = sum(Dose_info(:,1));
[T,Cv]=ode15s(@fun_DSB,[0 TotalDuration] ,c0 ); %time span
Cv_DSB = Cv(:,3);

    function dC = fun_DSB(time_ode,c)
        Cv_dNTP = [Cv_dNTP; c(1)];
        T_dNTP = [T_dNTP;  time_ode  ];
        [T_dNTP_unique, ia_dNTP, ic] = unique(T_dNTP ,'sorted');
        Cv_dNTP_2 =  Cv_dNTP(ia_dNTP);
        dcdt = zeros(3,1);
        

        
        PercentagefreeTS_q = interp1(T_UptoFreeTS , PercentagefreeTS , time_ode, 'PCHIP');
        
        theta_dNTP = theta(1:12);
        k1  =  theta_dNTP(1);
        k2 = theta_dNTP(2);
        k3 = theta_dNTP(3);
        k4 = theta_dNTP(4);
        k5 = theta_dNTP(5);
        k6 = theta_dNTP(6);
        k7 = theta_dNTP(7);
        k8  = theta_dNTP(8);
        k9 =  theta_dNTP(9);
        k_B =  theta_dNTP(10);
        k_A = theta_dNTP(11);
        k10  = theta_dNTP(12);
        gamma_dNTP = 0.15;
        
        
        TS_p_2 = (1-PercentagefreeTS_q).^gamma_dNTP;
        dcdt(1) =    k1*TS_p_2/( k2^gamma_dNTP +  TS_p_2)   - k3.*c(1)/(1+k4.*c(1)) - k5.*c(2)*c(1) /( (1+k6.*c(2))*(1+k7^gamma_dNTP .* TS_p_2) *( 1+k_B.*c(1) )   );
        %U
        dcdt(2) = k10*TS_p_2/(k8^gamma_dNTP+ TS_p_2)    - k9.*c(2)/(1+k_A.*c(2) )            ;

        Day2min = 24*60;
        theta_DSB = theta(13:end); 
        lag  = theta_DSB(1)*Day2min ; 
        V_dNTP = theta_DSB(2);
        K_dNTP =  theta_DSB(3);
        V_HRR =  theta_DSB(4);% theta_DSB(3);
        K_HRR = theta_DSB(5) ;%theta_DSB(4);
        k_i =  theta_DSB(6);
        k_0 =   theta_DSB(7); %theta(6);
        gamma_1  = 0.6;
        if time_ode <= lag
            dcdt(3)  = k_0 - ...
                V_HRR.*c(3) /(  K_HRR+c(3) );
        elseif time_ode > lag
             clag1 = interp1(T_dNTP_unique , Cv_dNTP_2  , time_ode-lag, 'PCHIP');
             dcdt(3)  = k_0.* (1+ V_dNTP.* clag1^gamma_1 / (  K_dNTP^gamma_1 + clag1^gamma_1  ))- ...
                V_HRR.*c(3)   /(  K_HRR+c(3) )/(k_i ^gamma_1* clag1^gamma_1  + 1  ) ;
        end
        dC=dcdt;
    end

end 

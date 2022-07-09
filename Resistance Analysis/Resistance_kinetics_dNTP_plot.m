function [T,Cv]=Resistance_kinetics_dNTP_plot(theta_cellular, theta_dNTP, theta_variation,Dose_info)
 
theta = theta_dNTP;
c0=[
    0
    0
    ];%initial condition
theta_resistance = theta_variation;
[T_TimeCourse_final,Cv_UPtoFreeTS_Timecourse_2]=Resistance_kinetics_uptoFreeTS(theta_cellular, theta_variation,Dose_info); %time span
Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:}; 
PercentagefreeTS = Cv_UPtoFreeTS_Timecourse_2(:,end);
[T,Cv]=ode15s(@DifEq,[0 Duration],c0);
Cv = real(Cv);
    function dC=DifEq(t_ode,c)
        PercentagefreeTS_q = interp1(T_TimeCourse_final,  PercentagefreeTS , t_ode, 'PCHIP'); 
        k1  =  theta(1);
        k2 = theta(2);
        k3 = theta(3);
        k4 = theta(4);
        k5 = theta(5);
        k6 = theta(6);
        k7 = theta(7);
        k8  = theta(8);
        k9 =  theta(9);
        k_B =  theta(10);
        k_A = theta(11);
        k10  = theta(12);
        % k10_min  = theta(13);
        %k1_min = theta(14);
        gamma_dNTP = 0.15;
        theta_resistance_cell = num2cell(theta_resistance);
       [ V_max_54,K_m_54,G_0, k_cat,k_08,k_95,k_59 ,k_09, k5,lag_DSB,alpha, V_HRR, K_HRR,  lag_tumor,E_max_damage, EC_50 ]  = theta_resistance_cell{:};
       

        
        dcdt = zeros(2,1);
        
        
        %TS_p  =   PercentagefreeTS_q.^gamma_dNTP  ;
        TS_p_2 = (1-PercentagefreeTS_q).^gamma_dNTP;
 
        %TS_p = abs(   TSf_0 /TS_0  -   c(9)/ (TS_0*( alpha + (1-alpha).*exp(-k_d.*t/60))  )   ) ./ TSf_0 /TS_0 ;
        %Nuc
        dcdt(1) =    k1*TS_p_2./( k2^gamma_dNTP +  TS_p_2)   - k3.*c(1)./(1+k4.*c(1)) - k5.*c(2)*c(1) ./( (1+k6.*c(2))*(1+k7^gamma_dNTP .* TS_p_2) *( 1+k_B.*c(1) )   );
        %U
        dcdt(2) = k10*TS_p_2./(k8^gamma_dNTP+ TS_p_2)    - k9.*c(2)./(1+k_A.*c(2) )            ;
        
        dC=dcdt;
    end

end
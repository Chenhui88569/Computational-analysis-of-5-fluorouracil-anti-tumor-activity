function [T,Cv_DSB] = RepeatedDoses_kinetics_dNTP_DSB(theta,theta_3to10, Model_Index,Control_flag)


AllModel_info = CollectDataIntoNestedStructure;
Dose_info = AllModel_info(Model_Index).Dose_regime;

Cv_dNTP  = [];
T_dNTP = [];
Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:}; 

c0=[
    0
    0
    42.755 % DSB
    ];%initial condition

%theta_3to10 = [0.969000000000000;540.079000000000;171.262000000000;537.075000000000;6.36200000000000;126.984000000000;25.1970000000000;2319.83400000000;0.102000000000000;0.490000000000000;0.0820000000000000;0.00600000000000000;0.0190000000000000;0.872000000000000;5.73400000000000;0.813000000000000;0.933000000000000;0.597000000000000;0.846000000000000;4.56300000000000;0.0330000000000000;2.85500000000000;0.108000000000000;21.8400000000000;0.168000000000000;17.4080000000000;0.0300000000000000];
%theta_3to10 =  [0.926830000000000;542.591670000000;172.642460000000;540.339100000000;6.96935000000000;133.508270000000;27.9590800000000;2324.40653000000;0.0708600000000000;0.503870000000000;0.0908900000000000;0.221610000000000;0.0186700000000000;0.851520000000000;5.51498000000000;0.958570000000000;0.755580000000000;0.536800000000000;0.284020000000000;4.34261000000000;0.0299300000000000;2.86930000000000;0.124640000000000;21.9376000000000;0.165060000000000;17.4138200000000;0.0269200000000000];
[T_UptoFreeTS , C_sim_UpToFreeTS ] = RepeatedDoses_kinetics_Cellular3to10(theta_3to10, Model_Index, Control_flag);
PercentagefreeTS = C_sim_UpToFreeTS(:,end);

[T,Cv]=ode15s(@fun_DSB,[0 Duration] ,c0 ); %time span
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

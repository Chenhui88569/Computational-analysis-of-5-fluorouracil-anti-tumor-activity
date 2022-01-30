function [T,Cv] = kinetics_DSB_plot(theta,Dose_info)

c0=[
    42.755 % DSB
    ];%initial condition
 
Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:}; 
theta_dNTP = [0.21351;32.50090;3.71407;73.34474;3.77027;0.30262;3.33454;1.35311;0.27614;35.66968;6.04530;0.07853];
[T_uptodNTP,Cv_dNTP]= kinetics_dNTP_plot(theta_dNTP,Dose_info);
dNTP_TimeCourse = Cv_dNTP(:,1);
[T,Cv]=ode15s(@fun_DSB,[0 Duration], c0 ); %time span
Cv = real(Cv);


    function dC = fun_DSB(time_ode,c)
        dcdt = zeros(1,1);
        
        lag  = theta(1) ; 
        V_dNTP = theta(2);
        K_dNTP =  theta(3);
        V_HRR =  theta(4);% theta(3);
        K_HRR = theta(5) ;%theta(4);
        k_i =  theta(6);
        k_0 =   theta(7); %theta(6);
        gamma_1  =0.6;
       
        if time_ode <= lag
            dcdt(1)  = k_0 - ...
                V_HRR.*c(1) /(  K_HRR+c(1) );
        elseif time_ode > lag
            clag1 = interp1(T_uptodNTP , dNTP_TimeCourse  , time_ode-lag, 'PCHIP');
             dcdt(1)  = k_0.* (1+ V_dNTP.* clag1^gamma_1 / (  K_dNTP^gamma_1 + clag1^gamma_1  ) )- ...
                V_HRR.*c(1)   /(  K_HRR+c(1) )/(k_i^gamma_1* clag1^gamma_1  + 1  ) ;
        end
        dC=dcdt;
    end

end 

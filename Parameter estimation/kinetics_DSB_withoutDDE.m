function Ysim = kinetics_DSB_withoutDDE(theta)
plasma_UnitConversion = 10^6/130.077; 

c0=[
    42.755 % DSB
    ];%initial condition

c_DSB = [42.755
74.758
74.276
113.8229%121.733
226.8139
189.0788%163.286
132.248
]; 

t_DSB_hour = [ 0
 8
24
32
48
64
72];
%t_DSB_hour = t_DSB_hour- 24*ones(5,1);
t_DSB = t_DSB_hour*60;

theta_dNTP = [0.21351;32.50090;3.71407;73.34474;3.77027;0.30262;3.33454;1.35311;0.27614;35.66968;6.04530;0.07853];
Dose_info = [ t_DSB(end)  0   0   93*plasma_UnitConversion  ];
[T_uptodNTP,Cv_dNTP]= kinetics_dNTP_plot(theta_dNTP,Dose_info);
dNTP_TimeCourse = Cv_dNTP(:,1);

[T,Cv]=ode15s(@fun_DSB, t_DSB ,c0 ); %time span
weight = ones( size(t_DSB_hour ));
weight(ismember(t_DSB_hour, [32  64])) = 0.00001;
%Ysim   = Cv(:,end)./c_DSB;
c_data_DSB_temp1 =c_DSB;
c_data_DSB_temp2= c_DSB ;
c_data_DSB_temp1 (c_data_DSB_temp1  ==0)=1;
c_output_DSB_normed = Cv ./c_data_DSB_temp1   ;
c_data_DSB_temp2( c_data_DSB_temp2~=0 )  =1;
Ysim = [c_output_DSB_normed- c_data_DSB_temp2].*weight; 

    function dC = fun_DSB(time_ode,c)
        dcdt = zeros(1,1);
      
        lag  = theta(1) ; 
        V_dNTP = theta(2);
        K_dNTP =  theta(3);
        V_HRR =  theta(4);% theta(3);
        K_HRR = theta(5) ;%theta(4);
        k_i =  theta(6);
        k_0 =   theta(7); %theta(6);
        gamma_1  = 0.6;
        if time_ode <= lag
            dcdt(1)  = k_0 - ...
                V_HRR.*c(1) /(  K_HRR+c(1) );
        elseif time_ode > lag
            clag1 = interp1(T_uptodNTP , dNTP_TimeCourse  , time_ode-lag, 'PCHIP');
             dcdt(1)  = k_0.* (1+ V_dNTP.* clag1^gamma_1 / (  K_dNTP^gamma_1 + clag1^gamma_1  ))- ...
                V_HRR.*c(1)   /(  K_HRR+c(1) )/(k_i ^gamma_1* clag1^gamma_1  + 1  ) ;
        end
        dC=dcdt;
    end

end 

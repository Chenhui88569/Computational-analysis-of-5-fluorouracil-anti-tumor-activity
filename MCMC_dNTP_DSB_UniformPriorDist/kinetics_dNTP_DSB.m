function Ysim = kinetics_dNTP_DSB(theta)

Cv_dNTP  = [];
T_dNTP = [];

c0=[
    0
    0
    42.755 % DSB
    ];%initial condition

c_DSB = [42.755
74.758
74.276
%113.8229%121.733
226.8139
%189.0788%163.286
132.248
]; 

t_DSB_hour = [ 0
 8
24
%32
48
%64
72];
%t_DSB_hour = t_DSB_hour- 24*ones(5,1);
t_DSB = t_DSB_hour*60;

t_dNTP_hour =[
       0 
0.5 
1
2
4
6
8
9
12 ];
t_dNTP = t_dNTP_hour.*60; % to minutes
Current_dNTP_M = 'alldNTP';
C_dNTP_measure_pertubation = dNTP_data_processing(Current_dNTP_M);    
Duration = [0 t_DSB(end)];
% C_dNTP_measure_pertubation = C_dNTP_measure_pertubation(non_add_point_idx );
%theta_3to10 = [0.980860000000000;538.450690000000;171.153100000000;536.641620000000;6.70454000000000;135.995300000000;27.1384400000000;2325.15421000000;0.0890600000000000;0.521000000000000;0.114470000000000;0.310230000000000;0.0199000000000000;0.839790000000000;5.53387000000000;0.959080000000000;0.807740000000000;0.554160000000000;0.395870000000000;4.35978000000000;0.0385200000000000;2.88205000000000;0.152570000000000;21.9351600000000;0.154750000000000;17.4304600000000;0.0321300000000000];
theta_3to10 = [0.84265;537.8213;180.865;530.619;6.9293;125.2409;24.4364;2323.3401;0.10407;0.49915;0.27889;0.14482;0.02128;0.82629;5.8241;0.81435;0.90936;0.56074;0.68855;4.4538;0.03442;2.8643;0.14933;21.9752;0.49934;17.4351;0.10468;200.0302;200.0388;5.0631;2.0197;0.29659;2.0088];
[T_UptoFreeTS , C_sim_UpToFreeTS ] = Plot_kinetics_Cellular3to10(theta_3to10,Duration  );
PercentagefreeTS = C_sim_UpToFreeTS(:,end);

%theta_dNTP = [0.21351;32.50090;3.71407;73.34474;3.77027;0.30262;3.33454;1.35311;0.27614;35.66968;6.04530;0.07853];
%[T_uptodNTP,Cv_dNTP]= Plot_kinetics_dNTP(theta_dNTP,Dose_info);
%dNTP_TimeCourse = Cv_dNTP(:,1);
theta_dNTP_DSB = theta;
[T,Cv]=ode15s(@fun_DSB, Duration ,c0 ); %time span
% weight = ones( size(t_DSB_hour ));
% weight(ismember(t_DSB_hour, [32  64])) = 0.00001;
%Ysim   = Cv(:,end)./c_DSB;
% c_data_DSB_temp1 =c_DSB;
% c_data_DSB_temp2= c_DSB ;
% c_data_DSB_temp1 (c_data_DSB_temp1  ==0)=1;
% c_output_DSB_normed = Cv ./c_data_DSB_temp1   ;
% c_data_DSB_temp2( c_data_DSB_temp2~=0 )  =1;
weight = ones( size(t_dNTP));
weight(ismember(t_dNTP_hour, [2,4,8,9])) = 0; 

Cv_dNTP_output = interp1(T , Cv(:,1)  ,t_dNTP, 'PCHIP');
Cv_DSB_output = interp1(T , Cv(:,3)  ,t_DSB, 'PCHIP');
Ysim = {  (Cv_dNTP_output -  C_dNTP_measure_pertubation(:,1)).*weight ;
   Cv_DSB_output  -  c_DSB };

    function dC = fun_DSB(time_ode,c)
        Cv_dNTP = [Cv_dNTP; c(1)];
        T_dNTP = [T_dNTP;  time_ode  ];
        [T_dNTP_unique, ia_dNTP, ic] = unique(T_dNTP ,'sorted');
        Cv_dNTP_2 =  Cv_dNTP(ia_dNTP);
        dcdt = zeros(3,1);
        

        
        PercentagefreeTS_q = interp1(T_UptoFreeTS , PercentagefreeTS , time_ode, 'PCHIP');
        
        theta_dNTP = theta_dNTP_DSB(1:12);
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
        theta_DSB = theta_dNTP_DSB(13:end); 
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

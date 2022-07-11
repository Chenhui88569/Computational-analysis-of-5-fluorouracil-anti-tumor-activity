function Ysim  =kinetics_dNTP(theta)
%para_num = 2 
plasma_UnitConversion = 10^6/130.077; 
day2min = 24*60;

c0=[
    0
    0
    ];%initial condition
t_dNTP_hour =[
       0 
0.5 
1
2
4
6
8
9
12 
  ];
% add_point = [2,4,8,9];
% non_add_point_idx = find(~ismember(t_dNTP_hour,add_point));
% t_dNTP_hour = t_dNTP_hour(non_add_point_idx);
t_dNTP = t_dNTP_hour.*60;
Current_dNTP_M = 'alldNTP';
C_dNTP_measure_pertubation = dNTP_data_processing(Current_dNTP_M);    
% C_dNTP_measure_pertubation = C_dNTP_measure_pertubation(non_add_point_idx );

Dose_info = [t_dNTP(end)	0	0	 93*plasma_UnitConversion]; 
[T_uptofreeTS, PercentagefreeTS_all]=kinetics_uptoFreeTS( Dose_info); %time span
PercentagefreeTS = PercentagefreeTS_all(:,end);
[T,Cv]=ode15s(@DifEq,t_dNTP ,c0);
weight = ones( size(t_dNTP));
weight(ismember(t_dNTP_hour, [2,4,8,9])) = 0.00001; 

Ysim = [Cv - C_dNTP_measure_pertubation].*weight;
% SSR =  sum( Ysim(:,1).^2);
% ObjFunction = @(num_datapoints, SSR, num_para) num_datapoints*log(SSR/num_datapoints)+2*(num_para+1)*num_datapoints...
%     /(num_datapoints-num_para-2);
% num_datapoints = size(Ysim,1);
% num_para =12;
% obj_output =   ObjFunction(num_datapoints,SSR,num_para);

    function dC=DifEq(t_ode,c)
        PercentagefreeTS_q = interp1(T_uptofreeTS , PercentagefreeTS , t_ode, 'PCHIP');
        
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
        gamma_dNTP = 0.15;
        dcdt = zeros(2,1);
        
        
        %TS_p  =   PercentagefreeTS_q.^gamma_dNTP ;
        TS_p_2 = (1-PercentagefreeTS_q).^gamma_dNTP;
        %TS_p = abs(   TSf_0 /TS_0  -   c(9)/ (TS_0*( alpha + (1-alpha).*exp(-k_d.*t/60))  )   ) ./ TSf_0 /TS_0 ;
        %Nuc
        %dcdt(1) =   (   k1_max - (k1_max- k1_min).*TS_p/( k2 + TS_p)  ) - k3.*c(1)/(1+k4.*c(1)) - k5.*c(2)*c(1) /( (1+k6.*c(2))*(1+k7.*TS_p_2) *( 1+k_B.*c(1) )   );
        dcdt(1) =    k1*TS_p_2/( k2^gamma_dNTP +  TS_p_2)   - k3.*c(1)/(1+k4.*c(1)) - k5.*c(2)*c(1) /( (1+k6.*c(2))*(1+k7^gamma_dNTP .* TS_p_2) *( 1+k_B.*c(1) )   );
        %U
        dcdt(2) = k10*TS_p_2/(k8^gamma_dNTP+ TS_p_2)    - k9.*c(2)/(1+k_A.*c(2) )            ;
        
        dC=dcdt;
    end

end
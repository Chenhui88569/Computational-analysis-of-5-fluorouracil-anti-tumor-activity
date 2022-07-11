function [Tsim,Csim,C_diff] = kinetics_plasma(theta_plasma,t_plasma)
%C_data =  cell2mat(Dataset_C_g{1});
c_plasma =  [ 92.9936306
70.2506231
67.0926138
48.329509
26.1156387 
11.3868339 
5.57423745
1.7525814]; % µg/mL  *10^6/130.077
c0_plasma=[  c_plasma(1) ;0];%initial condition %447880.98
[Tsim,Csim]=ode15s(@DifEq,t_plasma,c0_plasma);
t_plasma  = [0;10;20;30;60;90;120;240];
Csim_est = interp1(Tsim,Csim(:,1),t_plasma,'pchip'  );
C_diff =  Csim_est - c_plasma;
%F = C ./c_plasma; not scaled

    function dC=DifEq(~,c)
    dcdt=zeros(2,1);
    %k_e = theta_plasma(1);
    V_max = theta_plasma(1);
    Q_21 = theta_plasma(2);
    %k_12 = theta_plasma(3);
    %K_31 = theta_plasma(4);
    V_1	= theta_plasma(3);
    V_2	= theta_plasma(4);
    K_m = theta_plasma(5);
    %A = c.*[V_1,V_2];


    %dcdt(1) = Q_21/V_2 .*A(2)  - Q_21/V_1 .*A(1)-  V_max.* c(1)./(K_m+c(1)); %activity concentraion
    dcdt(1) = Q_21/V_1*c(2)  - Q_21/V_1 .*c(1)-  V_max.* c(1)./(K_m+c(1)); %activity concentraion
    %dcdt(1) = Q_21/V_2 .*A(2)  - Q_21/V_1 .*A(1)-  k_e*A(1); %activity concentraion
    dcdt(2) = Q_21/V_2.* c(1) - Q_21/V_2.*c(2); %activity concentraion
    dC=dcdt;
    end
end
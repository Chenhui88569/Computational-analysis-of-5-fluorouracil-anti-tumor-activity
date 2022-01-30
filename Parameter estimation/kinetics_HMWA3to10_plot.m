function [T,C]=kinetics_HMWA3to10_plot(theta,T)
%para_num = 2

c0=[92.9936306*10^6/130.077
    0
    0
    0
    0
    0
    0%7
    3.878 %dUMP final
    0 
   
    ];%initial condition
Cv_RNA  = [];
T_RNA = [];
Cv_DNA  = [];
T_DNA = [];
count_genome = 1;
%lag_RNA = 1*24*60-2*60 ;
%lag_DNA   = 3*24*60-2*60;
[T,Cv]=ode15s(@DifEq,T,c0);
alpha  =2.0212;
k_d =  1.0345;
TS_0 = 0.0186;
TS_total_output =  TS_0*(alpha + (1-alpha)*exp(-k_d.*T)  );
TS_free_output = TS_total_output -Cv(:,9);
TS_free_output  = reshape(  TS_free_output,  size(T)   );
C=  [ Cv(:,[3,4,5,6,7,8,9 ])   TS_free_output     ];
C(C<=0) =0;

    function dC=DifEq(t_ode,c)
        Cv_RNA(count_genome) = c(6);
        T_RNA(count_genome) = t_ode;
        [T_RNA_unique, ia_RNA, ic] = unique(T_RNA,'sorted');
        Cv_RNA_2 =  Cv_RNA(ia_RNA);
        
        Cv_DNA(count_genome) = c(7);
        T_DNA(count_genome) = t_ode;
        [T_DNA_unique, ia_DNA, ic] = unique(T_DNA,'sorted');
        Cv_DNA_2 =  Cv_DNA(ia_DNA);
        
        dcdt=zeros(9,1);
        %         theta_PK = [
        %             5609.02950
        %             4.81777
        %             379.09489
        %             791.41256
        %             90400.63504
        %             ];
        theta_PK = [5608.78414
0.01271
217.17977
453.41709
90394.64855];
        V_max_1 = theta_PK(1);
        Q_21 = theta_PK(2);
        V_1 = theta_PK(3);
        V_2 = theta_PK(4);
        K_m_1 = theta_PK(5);
        
        theta_FNUC = [ 539.9471708	171.0849025	536.8456805	6.34933249	126.555982	25.17731647	2321.536075	0.963912688
            ];
        
        V_max_influx =  theta_FNUC (1);
        K_m_influx =  theta_FNUC (2);
        V_max_efflux = theta_FNUC (3);
        K_m_efflux = theta_FNUC (4);
        k_03 = theta_FNUC (5);
        V_max_54 = theta_FNUC (6);
        K_m_54 = theta_FNUC (7);
        K_31 =   theta_FNUC (8); %theta14-21,originally
        %k_05 =theta_FNUC(9);
        %
%         theta_DNA = [
%             0.00005
%             1.49947
%             4.66608
%             243.56446]; %2
%         k_65  =0.013095072 ; %theta1,   0.013095072
%         k_56 =0.048178833; %theta2 ,0.048178833
%         V_max_75 = theta_DNA(1);
%         K_m_75 = theta_DNA(2);
%         V_max_7 = theta_DNA(3);
%         K_m_7 = theta_DNA(4);
%         
        k_65  = 0.013095072; %theta1,   0.013095072
        k_56 =0.048178833; %theta2 ,0.048178833
        theta_DNA = [0.01737
            0.86902
            5.78621
            0.96151
            ] ;
        V_max_75 =theta_DNA(1);
        K_m_75 = theta_DNA(2);
        V_max_7 =theta_DNA(3);
        K_m_7 = theta_DNA(4);
       
        theta_delaying = [1408.59448
            5415.53178
            0.01737
            0.86902
            5.78621
            0.96151];
        theta_delaying  = theta(1:6);
        lag_RNA = theta_delaying(1);
        lag_DNA =theta_delaying(2); 
        k_06 = theta_delaying(3);
        k_07 =  theta_delaying(4); 
        gamma_RNA =  theta_delaying(5);
        gamma_DNA =  theta_delaying(6);
        theta_TS = theta(7:end);
%          V_max_7 = theta(7);
%         K_m_7 =  theta(8);
        
%         theta_TS = [3.15430
%             1.38500
%             26.24083
%             0.21903
%             0.35653
%             1.31874
%             0.95754];

        K_dUMP =  theta_TS(1);
        G_0= theta_TS(2);
        k_cat =  theta_TS(3); 
        k08 =  theta_TS(4);
        k_95 =  theta_TS (5);
        k_59 =  theta_TS(6); 
        k_09 =  theta_TS(7);
         
        alpha  =2.0212;
        k_d =  1.0345;
        TS_0 = 0.0186;
        
      
        TS_total =  TS_0*(alpha + (1-alpha)*exp(-k_d.*t_ode)  );
        TS_free = TS_total - c(9);
        TS_dUMP =  c(8) *TS_free/( TS_free + K_dUMP) ;
        dcdt(1) = Q_21/V_2 .*c(2)- Q_21/V_1 .*c(1) - V_max_1.*c(1) /(K_m_1+c(1)); %activity concentraion
        dcdt(2) = Q_21/V_1.*c(1) - Q_21/V_2.*c(2);%activity concentraion
        dcdt(3) = K_31.*c(1) +V_max_efflux.*c(4) /(K_m_efflux+c(4))- V_max_influx.*c(3) /(K_m_influx+c(3))-...
            k_03.*c(3); % amount of activity
        dcdt(4) = -V_max_efflux.*c(4) /(K_m_efflux+c(4)) + V_max_influx.*c(3) /(K_m_influx+c(3))- ...
            V_max_54.*c(4) /(K_m_54+c(4));%amount of activity
        if  t_ode <= lag_RNA
            dcdt(5) =  V_max_54.*c(4)/(K_m_54+c(4) ) - k_65.*c(5) -V_max_75.*c(5)   /(K_m_75+c(5) ) ...
                - k_95.*(TS_total - TS_dUMP).*c(5) +k_56.*c(6)+k_59.*c(9) ;  %amount of activity
        elseif t_ode >= lag_RNA
            c_RNA_q = interp1(T_RNA_unique  ,  Cv_RNA_2  , t_ode-lag_RNA,  'PCHIP');
            dcdt(5) =  V_max_54.*c(4)/(K_m_54+c(4) ) - k_65.*c(5) -V_max_75.*c(5)   /(K_m_75+c(5) ) ...
                -k_95.*(TS_total - TS_dUMP).*c(5) + k_56.*c(6) /(1+ k_06.*c_RNA_q)^gamma_RNA+k_59.*c(9) ;  %amount of activity
        end
        if  t_ode <= lag_RNA
            dcdt(6) = k_65.*c(5) - k_56.*c(6);
       elseif  t_ode >  lag_RNA
          c_RNA_q = interp1(T_RNA_unique  ,  Cv_RNA_2  , t_ode-lag_RNA,  'PCHIP');
          % dcdt(6) = k_65.*c(5)- k_56.*c(6) /(1+ k_06.*c_RNA_q); %RNA
          dcdt(6) = k_65.*c(5)- k_56.*c(6) /(1+ k_06.*c_RNA_q)^gamma_RNA; %RNA
        end
        if  t_ode <= lag_DNA
            dcdt(7) = V_max_75.*c(5)  /(K_m_75+c(5) ) - V_max_7.*c(7)  /(K_m_7+c(7) );%DNA
            %                dcdt(7) = V_max_75.*c(5)  /(K_m_75+c(5) ) - k_07*c(7);%DNA
        elseif  t_ode >  lag_DNA
            c_DNA_q = interp1(T_DNA_unique  ,  Cv_DNA_2  , t_ode-lag_DNA,  'PCHIP');
            dcdt(7) = V_max_75.*c(5)  /(K_m_75+c(5) ) -V_max_7.*c(7)  /(K_m_7+c(7) ) /(1+ k_07*c_DNA_q )^gamma_DNA  ;%DNA
            %            dcdt(7) = V_max_75.*c(5)  /(K_m_75+c(5) ) - k_07*c(7) /(1+ k_07_lag*c_DNA_q )^gamma_DNA  ;%DNA
        end
       
       %total dUMP
        dcdt(8) =  G_0 - k_cat.*(  c(8)-c(8)/(1+TS_free /K_dUMP) ) - k08.*c(8)  ;
        %FdUMP-TS complex
        dcdt(9) = k_95.*(TS_total - TS_dUMP).*c(5) - k_59.*c(9)  -  k_09.*c(9) ;
        
        count_genome = count_genome+1;
        dC=dcdt;
    end


end
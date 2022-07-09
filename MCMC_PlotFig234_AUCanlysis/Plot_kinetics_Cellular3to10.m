function  [Tsim,Csim] = Plot_kinetics_Cellular3to10( theta,T_span)

t_genome = [ 0; 2;8; 24;48;72;168;240  ]*60;%h
c_RNA = [0
 9.61727039
6.487633256
9.155480343
2.365849466
2.880513238
0.338723071
0.041839875
] *10^3* 2*10^-3;
c_DNA = [0
   15.52683
13.606861
15.756227
9.8544088
14.428981
4.3345366
2.1127143  ]*0.2*10^-3; %
c_genome = [c_RNA' c_DNA']; %row vector

data_dUMP = [0  3.878
3 5.770
12 4.974
24 5.552
48 3.528] .*[60*ones(5,1)  ones(5,1)]; % time  value
data_TS_complete  =[   0         0    0.0186
    2   0.0355         0
    6    0.0307         0
   24   0.0272    0.0135
   32    0.0217    0.0161
   48    0.0180    0.0228].*[60*ones(6,1)  ones(6,1) ones(6,1)]; %time(h) , TSFdUMP_complex , free TS

 
c_TS_complete_vec =  [  transpose(data_TS_complete(:,2) )   transpose(data_TS_complete(:,3) ) ];%row vector

c_Intra5FU_FNUC= [0     0 
159.9649225	29.33645133
133.3177433	50.21923414
110.8097048	42.81788951
39.85969015	49.19029524
86.67641041	43.56620871
62.37942122	75.70885706
11.87956738	67.8631979
46.23209588	63.11604794
6.243788366	73.77959661
8.477053493	61.79479684
];
FNUC = c_Intra5FU_FNUC(:,2 );
Intracellular_5FU =  c_Intra5FU_FNUC(:,1 );
t_FNUC  = [0;15;30;40;50;60;70;80;90;105;120];% min, column vector
c_data = [Intracellular_5FU; FNUC ;c_genome'; c_TS_complete_vec'  ; data_dUMP(:,2)  ];% column vector
c0=[92.9936306
    0
    0
    0
    0
    0%RNA
    0% DNA
    data_dUMP(1,2) %dUMP final
    0%TS_FdUMP complex
 
      ];%initial condition 
  
Cv_RNA  = [];
T_RNA = [];
Cv_DNA  = [];
T_DNA = [];
count_genome = 1;
%lag_RNA = 1*24*60 -2*60;
%lag_DNA   = 3*24*60-2*60;  
[Tsim,Csim]=ode15s(@DifEq, T_span,c0);

alpha  =2.0212;
k_d =  1.0345;
TS_0 = 0.0186;

TS_total_output =  TS_0*(alpha + (1-alpha)*exp(-k_d.*Tsim)  );
TS_free_output = TS_total_output -  Csim(:,9) ;
Cv_FreeTS_Timecourse= TS_free_output./TS_total_output ; %add the free TS percentage to the last column
Csim(:,end+1) = Cv_FreeTS_Timecourse;
Csim(Csim<=0) =0;

% c_Intra5FU_output = interp1( T, Cv(:,4),  t_FNUC ,'pchip'  );
% c_FNUC_output = interp1( T, Cv(:,5),  t_FNUC ,'pchip'  );
% c_RNA_output = interp1( T, Cv(:,6),  t_genome ,'pchip'  );
% c_DNA_output = interp1( T, Cv(:,7),  t_genome ,'pchip'  );
% c_TScomplex_output = interp1( T, Cv(:,9),  data_TS_complete(:,1)  ,'pchip'  );
% %c_freeTS_output = interp1( T, Cv(:,9), data_TS_complete(:,1) ,'pchip'  );
% c_dUMP_output= interp1( T, Cv(:,8), data_dUMP(:,1)  ,'pchip'  );
% 
% c_output =   { c_Intra5FU_output- Intracellular_5FU;
%     c_FNUC_output - FNUC;
%     c_RNA_output-  c_RNA ;  
%     c_DNA_output-  c_DNA ; 
%     c_TScomplex_output  - data_TS_complete(:,2);
%        c_dUMP_output - data_dUMP(:,2)          };
%  

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
        theta_PK =   [329.539000000000,320.071000000000,16823.4720000000,2999.82100000000,17334.9770000000];
        
    
        V_max_1 = theta_PK(1);
        Q_21 = theta_PK(2);
        V_1 = theta_PK(3);
        V_2 = theta_PK(4);
        K_m_1 = theta_PK(5);
        
%         theta_FNUC = [ 539.9471708	171.0849025	536.8456805	6.34933249	126.555982	25.17731647	2321.536075	0.963912688
%             ];
        theta_FNUC = theta(1:8);
        K_31 =   theta_FNUC(1); %theta14-21,originally
        V_max_influx =  theta_FNUC(2);
        K_m_influx =  theta_FNUC(3);
        V_max_efflux = theta_FNUC(4);
        K_m_efflux = theta_FNUC(5);
        k_03 = theta_FNUC(6);
        V_max_54 = theta_FNUC(7);
        K_m_54 = theta_FNUC(8);
        
        %k_05 =theta_FNUC(9);
        %

        
        k_65  = theta(9); %0.013095072;  %RNA
        k_56 =  theta(10); %0.048178833;  %RNA
        k_06 = theta(11);
        gamma_RNA =  theta(12);
%         theta_DNA = [0.01737
%             0.86902
%             5.78621
%             0.96151] ;
      
        V_max_75 =  theta (13);
        K_m_75 =   theta(14);
        V_max_7 =  theta(15);
        K_m_7 =   theta (16);
        k_07 =  theta(17);
        gamma_DNA =  theta(18);
        lag_RNA =  theta(19)*24*60;
        lag_DNA = theta(20)*24*60;
  
        
        theta_TS = theta(21:end)  ;
%          V_max_7 = theta(7);
%         K_m_7 =  theta(8);
        
%         theta_TS = [3.15430
%             1.38500
%             26.24083
%             0.21903
%             0.35653
%             1.31874
%             0.95754];

        k_95 =  theta_TS (1);
        k_59 =  theta_TS(2); 
        k_09 =  theta_TS(3);
        K_dUMP =  theta_TS(4);
        G_0= theta_TS(5);
        k_cat =  theta_TS(6); 
        k08 =  theta_TS(7);
        
        
        alpha  =2.0212;
        k_d =  1.0345;
        TS_0 = 0.0186;
        
        TS_total =  TS_0*(alpha + (1-alpha)*exp(-k_d.*t_ode)  );
        TS_free = TS_total - c(9); 
        TS_dUMP = c(8) *TS_free/( TS_free + K_dUMP) ;
        dcdt(1) = Q_21/V_1*c(2)  - Q_21/V_1 .*c(1)-  V_max_1.* c(1)./(K_m_1 +c(1)); %activity concentraion
        %dcdt(1) = Q_21/V_2 .*A(2)  - Q_21/V_1 .*A(1)-  k_e*A(1); %activity concentraion
        dcdt(2) = Q_21/V_2.* c(1) - Q_21/V_2.*c(2); %activity concentraion
        dcdt(3) = K_31.*c(1)*10^6/130.077 +V_max_efflux.*c(4) /(K_m_efflux+c(4))- V_max_influx.*c(3) /(K_m_influx+c(3))-...
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
            dcdt(7) = V_max_75.*c(5)  /(K_m_75+c(5) ) -V_max_7.*c(7)  /(K_m_7+c(7) ) /(1+ k_07*c_DNA_q  )^gamma_DNA  ;%DNA
            %            dcdt(7) = V_max_75.*c(5)  /(K_m_75+c(5) ) - k_07*c(7) /(1+ k_07_lag*c_DNA_q )^gamma_DNA  ;%DNA
        end
         %dcdt(8) = k_85.*TS_total.*c(5) - k_58.*c(8)  -  k_08.*c(8) ;
        %total dUMP
        dcdt(8) =  G_0 - k_cat.*(  c(8)-c(8)/(1+TS_free /K_dUMP) ) - k08.*c(8)  ;
        %FdUMP-TS complex
        dcdt(9) = k_95.*(TS_total - TS_dUMP).*c(5) - k_59.*c(9)  -  k_09.*c(9) ;
       
%         dcdt(9) =  (  dcdt(11)  - dcdt(8) - c(9).* dcdt(10)/(     K_dUMP.*(c(9)/K_dUMP +1) ) )...
%             /(   c(10)/(c(9)+K_dUMP)  - c(10).*c(9)/( K_dUMP.^2.* (c(9)/K_dUMP +1).^2)  +1 )  ; 

      
        count_genome = count_genome+1;
        dC=dcdt;
    end


end
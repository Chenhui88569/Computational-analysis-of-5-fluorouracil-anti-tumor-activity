function  [T_TimeCourse_final,Cv_UPtoFreeTS_Timecourse_2] = RepeatedDoses_kinetics_Cellular3to10(theta,Model_Index,Control_flag)
warning('off','all')
AllModel_info = CollectDataIntoNestedStructure;
Dose_info = AllModel_info(Model_Index).Dose_regime;
if Control_flag == 1
    Dose_info(end) = 0;
end
% multiple doses 
Dose_info_cell = num2cell(Dose_info);
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:};
c0=[Dose
    0
    0
    0
    0
    0
    0%7
    3.878 %dUMP final
    0   
    ];%initial conditionCv_RNA  = [];
Cv_complex_Timecourse  = [];
Cv_complex_Timecourse(1,:) = c0' ;
Tv_TimeCourse =[];
Tv_TimeCourse(1) = 0; 

if Interval_TotalNum> 0
    for i =1:1: Interval_TotalNum
        %Initiate the time points 
        T_RNA = [];
        Cv_DNA  = [];
        Cv_RNA = [];
        T_DNA = [];
        count_genome = 1;
        %options = odeset('RelTol',1e-4,'AbsTol',1e-6);
        [Tv,Cv]=ode15s(@DifEq,0:100:DoseFrequency,c0 );
        Tv= Tv';
        Cv(Cv<=0) =0;
        c0(1) = Dose+ Cv(end,1);
        c0(2:end) = Cv(end,2:end);
        Tv_TimeCourse( end+1 : end+size(Tv,2)  )  = Tv + (i-1)*DoseFrequency;
        Cv_complex_Timecourse( end+1: end+size(Tv,2) ,:  ) =  Cv;
    end
end
%circument the situation where Duration=Interval_TotalNum*DoseFrequency
if Duration > Interval_TotalNum*DoseFrequency
    T_RNA = [];
    Cv_DNA  = [];
    Cv_RNA = [];
    T_DNA = [];
    count_genome = 1;
    t_interval_final = 0:0.1:Duration - Interval_TotalNum*DoseFrequency ;
   % options = odeset('RelTol',1e-4,'AbsTol',1e-6);
    %options = odeset('Stats','on');
    [Tv,Cv]=ode15s(@DifEq,t_interval_final,c0); %15s
%     else
%         [Tv,Cv]=ode15s(@DifEq,t_interval_final,c0); %15s
%     end
    Tv = Tv';
    Cv(Cv<=0) =0;
    Tv_TimeCourse( end+1 : end+size(Tv,2)  )  = Tv + Interval_TotalNum*DoseFrequency; %lump the time courses of interdose intervals
    Cv_complex_Timecourse( end+1: end+size(Tv,2),: ) =  Cv ;
end
 

[T_TimeCourse_final, ia_final, ~] = unique(Tv_TimeCourse,'sorted'); %horizontal
Cv_UPtoFreeTS_Timecourse_2 = Cv_complex_Timecourse( ia_final,:);

alpha  =2.0212;
k_d =  1.0345;
TS_0 = 0.0186;

TS_total_output =  TS_0*(alpha + (1-alpha)*exp(-k_d.*T_TimeCourse_final)  );
TS_free_output = TS_total_output - transpose(Cv_UPtoFreeTS_Timecourse_2(:,9));
Cv_FreeTS_Timecourse= TS_free_output./TS_total_output ; %add the free TS percentage to the last column
Cv_UPtoFreeTS_Timecourse_2(:,end+1) = Cv_FreeTS_Timecourse(:);
T_TimeCourse_final = T_TimeCourse_final';

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
        %theta_PK =   [329.539000000000,320.071000000000,16823.4720000000,2999.82100000000,17334.9770000000];
        
        theta_PK = theta(1:5);
        V_max_1 = theta_PK(1);
        Q_21 = theta_PK(2);
        V_1 = theta_PK(3);
        V_2 = theta_PK(4);
        K_m_1 = theta_PK(5);
        
%         theta_FNUC = [ 539.9471708	171.0849025	536.8456805	6.34933249	126.555982	25.17731647	2321.536075	0.963912688
%             ];
        theta_FNUC = theta(5+1:5+8);
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

        
        k_65  = theta(5+9); %0.013095072;  %RNA
        k_56 =  theta(5+10); %0.048178833;  %RNA
        k_06 = theta(5+11);
        gamma_RNA =  theta(5+12);
%         theta_DNA = [0.01737
%             0.86902
%             5.78621
%             0.96151] ;
      
        V_max_75 =  theta(5+13);
        K_m_75 =   theta(5+14);
        V_max_7 =  theta(5+15);
        K_m_7 =   theta(5+16);
        k_07 =  theta(5+17);
        gamma_DNA =  theta(5+18);
        lag_RNA =  theta(5+19);
        lag_DNA = theta(5+20);
  
        
        theta_TS = theta(5+21:end)  ;
%          V_max_7 = theta(7);
%         K_m_7 =  theta(8);
        
%         theta_TS = [3.15430
%             1.38500
%             26.24083
%             0.21903
%             0.35653
%             1.31874
%             0.95754];

        k_95 =  theta_TS(1);
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
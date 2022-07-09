function  Para_permutated = CreateParaPopulation(Anabolites_flag, Anabolites_TestedPara,...
                TSInhibition_flag ,TSInhibition_TestedPara,...
                dNTP_flag,DSB_flag, TGI_flag,DSB_TestedPara, TGI_TestedPara, ParaPopulation_size )
  
%CREATEPARAPOPULATION Summary of this function goes here
%   Detailed explanation goes here
%% assign parameter value to the submodels
load('Para.mat','Para_Info');
Para_value = Para_Info.Para_Value;
theta_cellular = Para_value(1:32);
%theta_FNUC = [25.17731647	2321.536075];
V_max_54 = theta_cellular(12);
K_m_54 = theta_cellular(13);
theta_TS = [theta_cellular(26:end);2.0212;1.0345;0.0186]; 
theta_TS_cell = num2cell(theta_TS);
[k_95, k_59 ,k_09, K_dUMP, G_0, k_cat, k_08, alpha, k_d, TS_0] = theta_TS_cell{:}; 

theta_dNTP = [ Para_value(33:33+11); 0.15 ];
theta_dNTP_cell = num2cell(theta_dNTP);
[k1, k2,k3,k4,k5,k6,k7,k8,k9, k_B,k_A,k10 ,gamma_dNTP ] = theta_dNTP_cell{:};

theta_DSB = Para_value(33+11+1:33+18);
theta_DSB_cell = num2cell(theta_DSB);
[lag_DSB , V_dNTP, K_dNTP, V_HRR, K_HRR, k_i,k_0] =  theta_DSB_cell{:}; 

theta_TGI_50_sensitive = [1339.31457867841;104.059993350524;0.841096217568630;73.9835830977462;0.000227131913701615;11.2628519108548;4.57756840800819e-06];
theta_TGI_cell = num2cell(theta_TGI_50_sensitive);
[lag_tumor, IC50,E_max_damage, EC_50 ,lambda_g , P_max, lambda_d ] = theta_TGI_cell{:}; 

ChangedPara = [];
if Anabolites_flag ==1
    ChangedPara = [ChangedPara  Anabolites_TestedPara  ];
    if find(Anabolites_TestedPara == 'V_max_54' )
       LowerBound_fold  = 2;
       ub  = V_max_54; 
       lb = 1/LowerBound_fold*V_max_54;
       Population_Vm54 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
       Population_Vm54= V_max_54*ones(ParaPopulation_size ,1 );
    end
    if find(Anabolites_TestedPara =='K_m_54')
        UpperBound_fold = 2;
        ub = UpperBound_fold* K_m_54;
        lb =  K_m_54;
        Population_Km54 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_Km54 =   K_m_54 *ones(ParaPopulation_size ,1 );
    end
else
     SplitOnCol = num2cell( ones(ParaPopulation_size,1 ) * [ V_max_54,  K_m_54  ],1);%1 :split on column
     [ Population_Vm54,Population_Km54] = SplitOnCol{:}; 
end

if dNTP_flag == 1 % k5 increse the rate of recovery  
    ChangedPara = [ChangedPara "k5"  ];
    UpperBound_fold = 10 ;
    lb =  k5 ;
    ub =UpperBound_fold*k5 ;
    Population_k5 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
else
    Population_k5 =  k5*ones(ParaPopulation_size ,1 ); %vector
end

if TSInhibition_flag == 1 %"G_0","cat","08"  "59","95","09"
    ChangedPara = [ChangedPara TSInhibition_TestedPara ];
    if find(TSInhibition_TestedPara == 'G_0' )
        UpperBound_degree = 10;
        lb =   G_0;
        ub = UpperBound_degree *G_0;
        Population_G0 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_G0 =  G_0*ones(ParaPopulation_size ,1 );
    end 
    if find(TSInhibition_TestedPara =='cat')
        LowerBound_fold = 10 ;
        lb =  1/LowerBound_fold*k_cat;
        ub = LowerBound_fold * k_cat;
        Population_cat =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_cat =  k_cat*ones(ParaPopulation_size ,1 );
    end
    if find(TSInhibition_TestedPara =='08') %dUMP accumulation
        LowerBound_fold = 10;
        ub =   k_08;
        lb = 1/LowerBound_fold *k_08;
        Population_08 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_08 = k_08*ones(ParaPopulation_size ,1 );
    end
    if find(TSInhibition_TestedPara =='59')%dissociation
        UpperBound_fold = 10 ;
        lb =  k_59;
        ub =UpperBound_fold*k_59;
        Population_59 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_59 =  k_59*ones(ParaPopulation_size ,1 );
    end
    if find(TSInhibition_TestedPara =='95') %association
        LowerBound_fold = 10;
        lb = 1/LowerBound_fold*k_95;
        ub = k_95;
        Population_95 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_95 =  k_95*ones(ParaPopulation_size ,1 );%retain value
    end
    if find(TSInhibition_TestedPara =='09') %degradation/elimination 
        UpperBound_fold = 10;
        ub = UpperBound_fold* k_09;
        lb =  k_09;
        Population_09 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_09 =  k_09*ones(ParaPopulation_size ,1 );%retain value
    end
    
    if find(TSInhibition_TestedPara == 'alpha' )
        UpperBound_degree = 10;
        lb =  alpha;
        ub = UpperBound_degree *alpha;
        Population_alpha =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_alpha =  alpha*ones(ParaPopulation_size ,1 );
    end
else
    %if entire TS mechanism is not of interest
   SplitOnCol = num2cell( ones(ParaPopulation_size,1 ) * [ G_0, k_cat,k_08,k_95,k_59 ,k_09,alpha   ] ,1);%1 :split on column
   [ Population_G0, Population_cat,Population_08 ,Population_95 , Population_59, Population_09,  Population_alpha]  = SplitOnCol{:};
end

if DSB_flag == 1 
    ChangedPara = [ChangedPara  DSB_TestedPara   ];
    if find( DSB_TestedPara == 'lag_DSB' )
         UpperBound_fold = 10 ;
         lb =lag_DSB   ;
         ub =UpperBound_fold * lag_DSB   ;
         Population_Tdsb =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_Tdsb =  lag_DSB  *ones(ParaPopulation_size ,1 ); %vector
    end
    if find( DSB_TestedPara == 'V_HRR' )
        UpperBound_fold = 10 ;
        lb = V_HRR   ;
        ub =UpperBound_fold * V_HRR    ;
        Population_VHR =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_VHR  =  V_HRR *ones(ParaPopulation_size ,1 ); %vector
    end
    
    if find( DSB_TestedPara == 'K_HRR' )
        LowerBound_fold = 10;
        lb = 1/LowerBound_fold * K_HRR;
        ub = K_HRR ;
        Population_KHR =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_KHR =   K_HRR *ones(ParaPopulation_size ,1 ); %vector
    end
   
else 
     SplitOnCol = num2cell( ones(ParaPopulation_size,1 ) * [ lag_DSB ,  V_HRR, K_HRR  ],1);%1 :split on column
     [Population_Tdsb ,Population_VHR, Population_KHR] = SplitOnCol{:};
end

if TGI_flag == 1 
    ChangedPara = [ChangedPara  TGI_TestedPara    ];
    if find(  TGI_TestedPara == 'lag_tumor_TGI' )
         UpperBound_fold = 10 ;
         lb = lag_tumor ;
         ub =UpperBound_fold* lag_tumor ;
         Population_TdTGI =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_TdTGI =  lag_tumor * ones(ParaPopulation_size ,1 ); %vector
    end
    if find(  TGI_TestedPara == 'E_max_damage' )
        LowerBound_fold = 10;
        lb = 1/LowerBound_fold * E_max_damage;
        ub = E_max_damage ;
        Population_Edamage =  lb + ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_Edamage   =  E_max_damage *ones(ParaPopulation_size ,1 ); %vector
    end
    
    if find(  TGI_TestedPara == 'EC_50_damage' )
        UpperBound_fold = 10 ;
        lb = EC_50  ;
        ub =UpperBound_fold * EC_50   ;
        Population_EC50 =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    );
    else
        Population_EC50 =   EC_50 * ones(ParaPopulation_size ,1 ); %vector
    end
else 
     SplitOnCol = num2cell( ones(ParaPopulation_size,1 ) * [ lag_tumor,E_max_damage, EC_50 ],1);%1 :split on column
     [Population_TdTGI , Population_Edamage , Population_EC50] = SplitOnCol{:};
end

Para_unpermutated = [ Population_Vm54   Population_Km54  Population_G0 Population_cat  Population_08   Population_95...
    Population_59 Population_09  Population_k5  Population_Tdsb , Population_alpha,  Population_VHR,  Population_KHR,...
    Population_TdTGI  Population_Edamage   Population_EC50];
Para_permutated = Para_unpermutated;
end


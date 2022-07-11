function Target_value = Target_fun_LE_TGI(theta,Model_Index)
% if we treat the measurement error as random variable. 
% and it's time-independent.
C_diff_cell =  kinetics_TumorVolume_Treated(theta,Model_Index);
num_dataset = size(C_diff_cell,1);
ThetaNoise = theta(end-num_dataset +1 :end);
Target_value  = 1 ;
%% specify the prior distribution for each parameter in the TGI model
% %% Treated for 100mg/kg
if Model_Index  == 1
    prior_T_d_Tv                = normpdf(theta(1),500,100);
    prior_IC50                       = wblpdf(theta(2),10,2);
    prior_E_max_damage   = unifpdf(theta(3),1e-8,1);
    prior_EC_50_damage      = wblpdf(theta(4),5.2,2);
    prior_g                         = unifpdf(theta(5),1e-8,1);%1e-3
    prior_Pmax                      = wblpdf(theta(6),5,2);
    prior_d                         = unifpdf(theta(7),1e-8,1); %1e-3
    % prior_nosie_treated     =  lognpdf(ThetaNoise(1),4,0.4);
    % prior_nosie_control     =    gampdf( ThetaNoise(2),4,3 );
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );
elseif Model_Index  == 2
    % Treated for 50mg/kg
    prior_T_d_Tv                = normpdf(theta(1),500,100);
    prior_IC50                       =  unifpdf(theta(2),1e-8,5); %1e-3
    prior_E_max_damage   = wblpdf(theta(3),22,2) ;
    prior_EC_50_damage      = wblpdf(theta(4),157,2);
    prior_g                         = unifpdf(theta(5),1e-8,1);%1e-3
    prior_Pmax                      =  unifpdf(theta(6),1,10); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(7),1e-8,1); %1e-3
    % prior_nosie_treated     =  lognpdf(ThetaNoise(1),4,0.4);
    % prior_nosie_control     =    gampdf( ThetaNoise(2),4,3 );
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );
    
elseif Model_Index  == 3
    % Treated for 20 mg/kg Twice a week for four weeks
    prior_T_d_Tv                = normpdf(theta(1),900,100);
    prior_IC50                       =  unifpdf(theta(2),1e-8,2); %enlarge from 1
    prior_E_max_damage   =  unifpdf(theta(3),1e-8,2);
    prior_EC_50_damage      = unifpdf(theta(4),1e-8,2);
    prior_g                         = unifpdf(theta(5),1e-8,1);%1e-3
    prior_Pmax                      =  unifpdf(theta(6),1,2); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(7),1e-8,1); %1e-3
    % prior_nosie_treated     =  lognpdf(ThetaNoise(1),4,0.4);
    % prior_nosie_control     =    gampdf( ThetaNoise(2),4,3 );
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );
    
elseif Model_Index  == 4
    prior_T_d_Tv = 1;
    prior_IC50 = 1;
    prior_E_max_damage = 1;
    prior_EC_50_damage = 1;
    prior_g                         = unifpdf(theta(1),1e-8,1);%1e-3
    prior_Pmax                      =  unifpdf(theta(2),1,5); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(3),1e-8,1); %1e-3
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );

elseif Model_Index  == 5
    % Treated for 20 mg/kg Twice a week for four weeks
    prior_T_d_Tv                = normpdf(theta(1),800,100);
    prior_IC50                       =  unifpdf(theta(2),1e-8,2); %enlarge from 1
    prior_E_max_damage   =  unifpdf(theta(3),1e-8,2);
    prior_EC_50_damage      = unifpdf(theta(4),1e-8,2);
    prior_g                         = unifpdf(theta(5),1e-8,1);%1e-3
    prior_Pmax                      =  unifpdf(theta(6),0.5,4); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(7),1e-8,1); %1e-3
    % prior_nosie_treated     =  lognpdf(ThetaNoise(1),4,0.4);
    % prior_nosie_control     =    gampdf( ThetaNoise(2),4,3 );
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );
elseif Model_Index  == 6
    prior_T_d_Tv = 1;
    prior_IC50 = 1;
    prior_E_max_damage = 1;
    prior_EC_50_damage = 1;
    prior_g                         = unifpdf(theta(1),1e-8,1);%1e-3
    prior_Pmax                      =  wblpdf(theta(2),13,2) ; %unifpdf(theta(2),8,16); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(3),1e-8,1); %1e-3
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );

 elseif Model_Index  == 7
    prior_T_d_Tv = 1;
    prior_IC50 = 1;
    prior_E_max_damage = 1;
    prior_EC_50_damage = 1;
    prior_g                         = unifpdf(theta(1),1e-8,1);%1e-3
    prior_Pmax                      =  unifpdf(theta(2),0.5,4) ; %unifpdf(theta(2),8,16); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(3),1e-8,1); %1e-3
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );
elseif Model_Index  == 8
    prior_T_d_Tv                = normpdf(theta(1),800,100);
    prior_IC50                       = normpdf(theta(2),50,100);
    prior_E_max_damage   = unifpdf(theta(3),1e-8,1);
    prior_EC_50_damage    = unifpdf(theta(4),1e-8,100);
    prior_g                         = unifpdf(theta(5),1e-8,1);%1e-3
    prior_Pmax                      =  wblpdf(theta(6),5,2) ; %unifpdf(theta(2),8,16); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(7),1e-8,1); %1e-3
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );
elseif Model_Index == 9
    prior_T_d_Tv                = normpdf(theta(1),900,100);
    prior_IC50                       = normpdf(theta(2),50,100);
    prior_E_max_damage   = unifpdf(theta(3),1e-8,5);
    prior_EC_50_damage    = unifpdf(theta(4),1e-8,100);
    prior_g                         = unifpdf(theta(5),1e-8,1);%1e-3
    prior_Pmax                      =  wblpdf(theta(6),5,2) ; %unifpdf(theta(2),8,16); %1e-3 %wblpdf(theta(6),3,2);
    prior_d                         = unifpdf(theta(7),1e-8,1); %1e-3
    prior_nosie_treated     =  gampdf( ThetaNoise(1), 1,1 );
    prior_nosie_control     =  gampdf( ThetaNoise(2), 1,1 );

end



prior_vec = [
    prior_T_d_Tv
    prior_IC50
    prior_E_max_damage
    prior_EC_50_damage
    prior_g
    prior_Pmax
    prior_d
    prior_nosie_treated
    prior_nosie_control
    ];
prior_val  = prod(prior_vec);

%% calculate likehood for each dataset
for i = 1:num_dataset
    temp_diff = C_diff_cell{i};
    num_Data = length( temp_diff);
    diff_norm = vecnorm(temp_diff,2).^2;
    Noise_Term = ThetaNoise(i); %variance
    C_vec = (1./ (2*pi*Noise_Term) ).^(num_Data./2); %vector % constant factor 2pi doesn't influence the result.
    Target_LE_temp = C_vec*...
    exp(- 1/(2*Noise_Term) *diff_norm)    ;  % eliminate the multiplicative constant
    Target_value =   Target_value*Target_LE_temp;
end

%% multiply likehood with prior
Target_value   = Target_value *prior_val;
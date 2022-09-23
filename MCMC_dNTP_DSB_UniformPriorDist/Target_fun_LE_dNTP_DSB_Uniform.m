function Target_value = Target_fun_LE_dNTP_DSB_Uniform(theta)
% if we treat the measurement error as random variable. 
% and it's time-independent.
C_diff_cell =  kinetics_dNTP_DSB(theta); 
num_dataset = size(C_diff_cell,1);
ThetaNoise = theta(end-num_dataset +1 :end);
Target_value  = 1 ;
%% specify the prior distribution for each parameter in the cellular mode
%% dNTP
% The variances of normal random variable used to 2.
prior_k_1 = unifpdf(theta(1),1e-3,1);
%prior_k_1 = gampdf(theta_dNTP(1),1,1);
prior_k_2 =  unifpdf(theta(2),0.1,100);  %normpdf(theta(2),32,10);% 
prior_k_3 = unifpdf(theta(3),0.01,10); %wblpdf(theta(3),4,2);
prior_k_4 =  unifpdf(theta(4),0.1,200);  %wblpdf(theta(4),73,2);
prior_k_5 = unifpdf(theta(5),0.01,10); %wblpdf(theta(5),4,2);
%prior_k_6 = gampdf(theta(6),1,1);
prior_k_6 = unifpdf(theta(6),1e-3,1);
prior_k_7 = unifpdf(theta(7),0.01,10); %wblpdf(theta(7),4,2);
prior_k_8 = unifpdf(theta(8),0.01,50); %wblpdf(theta(8),10,1); 
prior_k_9 = unifpdf(theta(9),1e-3,1);
prior_k_B =  unifpdf(theta(10),1,100);  %normpdf(theta(10),35,10);
prior_k_A =   unifpdf(theta(11),0.1,20);  %wblpdf(theta(11),7,2);
prior_k_10 = unifpdf( theta(12),1e-4,0.4) ;
prior_nosie_dNTP      =   unifpdf( ThetaNoise(1),0.001, 10) ; %gampdf( ThetaNoise(1), 1,1 );
%% DSB
prior_Td_DSB =  unifpdf( theta(13),0.1, 5) ; %wblpdf(theta(13),1.8,2);
prior_Vmax_dNTP =   unifpdf( theta(14),0.1, 50) ; %wblpdf(theta(14),10,2);
prior_Km_dNTP    = unifpdf( theta(15),1, 500) ; %normpdf(theta(15),152,50);
prior_Vmax_HR    =  unifpdf( theta(16),0.1, 50) ;    %wblpdf(theta(16),10,2);
prior_Km_HR       =  unifpdf( theta(17),1, 500) ;   %normpdf(theta(17),194,50);
prior_k_i           =  unifpdf( theta(18),1e-4,1) ; %0.5
prior_k_0           = unifpdf( theta(19),0.1,10) ;  %wblpdf(theta(19),3,2);
prior_nosie_DSB      =  unifpdf( ThetaNoise(2),5,500) ;% lognpdf(ThetaNoise(2),4,0.4);



prior_vec = [
    prior_k_1
    prior_k_2
    prior_k_3
    prior_k_4
    prior_k_5
    prior_k_6
    prior_k_7
    prior_k_8
    prior_k_9
    prior_k_B
    prior_k_A
    prior_k_10
    prior_Td_DSB
    prior_Vmax_dNTP
    prior_Km_dNTP
    prior_Vmax_HR
    prior_Km_HR
    prior_k_i
    prior_k_0
    prior_nosie_dNTP
    prior_nosie_DSB
    
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
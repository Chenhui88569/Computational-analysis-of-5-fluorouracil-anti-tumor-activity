function Target_value = Target_fun_LE_Cellular(theta)
% if we treat the measurement error as random variable. 
% and it's time-independent.
C_diff_cell = kinetics_Cellular3to10(theta);
num_dataset = size(C_diff_cell,1);
ThetaNoise = theta(end-num_dataset+1:end);
Target_value = 1;

%% specify the prior distribution for each parameter in the cellular mode

% %% %5-FU in interstitial fluid
% prior_Q31 = wblpdf(theta(1) ,1,2); 
% prior_Vinflux = normpdf(theta(2) , 5.4e2,2);
% prior_Kinflux = normpdf(theta(3) , 1.71e2,2);
% prior_Vefflux = normpdf(theta(4) , 5.37e2,2);
% prior_Kefflux = normpdf(theta(5) , 6.35,2);
% prior_k03 =  normpdf( theta(6), 1.27e2,2);%chi2pdf(theta(9:10),1);
% prior_V54 = normpdf(theta(7), 2.52e1,2);
% prior_K54 = normpdf(theta(8), 2.32e3,2);
% %% 5FU incorporation into DNA and RNA
% %prior_k_65 = gampdf(theta(9), 1,1 )  ;
% prior_k_65 = unifpdf( theta(9),1e-3,0.8) ;
% %prior_k_56 = gampdf(theta(10), 1,1 )  ;
% prior_k_56 = unifpdf(theta(10),1e-3,1.1 )  ;
% %prior_k_06 = gampdf(theta(11), 1,1 ) ;
% prior_k_06 = unifpdf( theta(11),1e-5,0.8) ;
% %prior_gamma_RNA = gampdf(theta(12), 1,1 ) ;
% prior_gamma_RNA = unifpdf(theta(12), 1e-5,0.8) ;
% prior_V_75 = gampdf(theta(13), 1,1);
% prior_K_75 = gampdf(theta(14), 1,1);
% prior_V_57 = wblpdf(theta(15), 5,2);
% prior_K_57 = wblpdf(theta(16), 1,2);
% prior_k07   = wblpdf(theta(17), 1,2);
% prior_gamma_DNA = wblpdf(theta(18), 1,2);
% prior_T_RNA = wblpdf(theta(19), 1,2);
% prior_T_DNA = wblpdf(theta(20), 5,2);
% prior_k95 = gampdf( theta(21), 1,1       );
% % prior_k95  = unifpdf( theta(21),1e-3,0.2) ;
% prior_k59 = wblpdf( theta(22), 5,2       );
% prior_k09 = gampdf( theta(23), 1,1     );
% 
% prior_KdUMP = wblpdf( theta(24), 17,5    );
% prior_G0 = gampdf( theta(25), 1,1    );
% prior_kcat = wblpdf( theta(26), 16,4  );
% %prior_k08 = gampdf( theta(27),1,1 ) ;
% prior_k08 = unifpdf( theta(27),1e-3,0.8) ;

%% %5-FU in interstitial fluid
prior_Q31 = wblpdf(theta(1) ,1,2); 
prior_Vinflux = normpdf(theta(2) , 5.4e2,50);
prior_Kinflux = normpdf(theta(3) , 1.71e2,20);
prior_Vefflux = normpdf(theta(4) , 5.37e2,50);
prior_Kefflux = normpdf(theta(5) , 6.35,3);
prior_k03 =  normpdf( theta(6), 1.27e2,20);%chi2pdf(theta(9:10),1);
prior_V54 = normpdf(theta(7), 2.52e1,10);
prior_K54 = normpdf(theta(8), 2.32e3,50);
%% 5FU incorporation into DNA and RNA
%prior_k_65 = gampdf(theta(9), 1,1 )  ;
prior_k_65 = unifpdf( theta(9),1e-3,0.8) ;
%prior_k_56 = gampdf(theta(10), 1,1 )  ;
prior_k_56 = unifpdf(theta(10),1e-3,1.1 )  ;
%prior_k_06 = gampdf(theta(11), 1,1 ) ;
prior_k_06 = unifpdf( theta(11),1e-5,0.8) ;
%prior_gamma_RNA = gampdf(theta(12), 1,1 ) ;
prior_gamma_RNA = unifpdf(theta(12), 1e-5,0.8) ;
prior_V_75 = gampdf(theta(13), 1,1);
prior_K_75 = gampdf(theta(14), 1,1);
prior_V_57 = wblpdf(theta(15), 5,2);
prior_K_57 = wblpdf(theta(16), 1,2);
prior_k07   = wblpdf(theta(17), 1,2);
prior_gamma_DNA = wblpdf(theta(18), 1,2);
prior_T_RNA = wblpdf(theta(19), 1,2);
prior_T_DNA = wblpdf(theta(20), 5,2);
prior_k95 = gampdf( theta(21), 1,1       );
% prior_k95  = unifpdf( theta(21),1e-3,0.2) ;
prior_k59 = wblpdf( theta(22), 5,2       );
prior_k09 = gampdf( theta(23), 1,1     );

prior_KdUMP = wblpdf( theta(24), 17,2    );
prior_G0 = gampdf( theta(25), 1,1    );
prior_kcat = wblpdf( theta(26), 18,2  );
%prior_k08 = gampdf( theta (27),1,1 ) ;
prior_k08 = unifpdf( theta(27),1e-3,0.8) ;

prior_noise_intra    = lognpdf(ThetaNoise(1),4,0.4);
prior_nosie_FNUC  =  lognpdf(ThetaNoise(2),4,0.4); 
prior_nosie_RNA   =   gampdf( ThetaNoise(3),4,1 );
prior_nosie_DNA   =  gampdf( ThetaNoise(4), 1,1 );
prior_nosie_TS      =  gampdf( ThetaNoise(5), 1,1 );
prior_nosie_dUMP =  lognpdf(ThetaNoise(6),1,0.5); 
prior_vec = [
    prior_Q31 
    prior_Vinflux
    prior_Kinflux 
    prior_Vefflux 
    prior_Kefflux
    prior_k03
    prior_V54
    prior_K54
    prior_k_65
    prior_k_56
    prior_k_06
    prior_gamma_RNA
    prior_V_75
    prior_K_75
    prior_V_57
    prior_K_57
    prior_k07
    prior_gamma_DNA
    prior_T_RNA
    prior_T_DNA
    prior_k95
    prior_k59
    prior_k09
    prior_KdUMP
    prior_G0
    prior_kcat
    prior_k08
    prior_noise_intra
    prior_nosie_FNUC
    prior_nosie_RNA
    prior_nosie_DNA
    prior_nosie_TS
    prior_nosie_dUMP
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
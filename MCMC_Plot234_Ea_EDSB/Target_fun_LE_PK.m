function Target_LE_PK = Target_fun_LE_PK(theta)
% if we treat the measurement error as random variable. 
% and it's time-independent.
t_plasma  = [0;10;20;30;60;90;120;240];
[~,~,C_diff] = kinetics_plasma(theta,t_plasma);
diff_norm = vecnorm(C_diff,2).^2;
Noise_Term = theta(end);
num_Data = length(C_diff);
C_vec = (1./Noise_Term).^(num_Data./2);
%rand_noise = 0.2;  % the measurement noise is kept as constant
Target_LE_PK = prod(C_vec)*...
    exp(- 1/(2*Noise_Term) *diff_norm)    ;  % eliminate the multiplicative constant
   
end
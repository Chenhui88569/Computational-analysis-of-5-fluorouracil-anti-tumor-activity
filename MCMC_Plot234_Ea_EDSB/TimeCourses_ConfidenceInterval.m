function  CI = TimeCourses_ConfidenceInterval(...
   KineticsPopulation,  Kinetics_Mean, ParaPopulation_size    ) %CI: confidence interval
%  the dimension of the KineticsPopulation is N*Np, N is the number of
%  parameters sets, Np is the the number of time points .
%  std(KineticsPopulation ) retruns the standard deviation of each time
%  point.
dbstop if error
% SEM_TS = std(KineticsPopulation )/sqrt(ParaPopulation_size  ); %row vector 
% conf = 0.95;
% alpha  = 1- conf;
% pLo = alpha  /2;
% pUp = 1- pLo ;
% ts = tinv( [pLo  pUp ], ParaPopulation_size-1 ); %95% confidence level 
% %c_mean = [ transpose(mean( KineticsPopulation))  transpose(mean( KineticsPopulation))];%the mean of each column
% CI = Kinetics_Mean + SEM_TS'*ts;  %% time course mean +
CI = quantile(KineticsPopulation,[0.15,0.85]);
CI = CI';
CI(CI<0) = 0;
%C_mean = c_mean(:,1); %vertical
%find the overall max ;the maximum minimum
%zero crossing can not be applied to all the situations occurring in the
%model
%compute the distance matrix  
% D_vec = vecnorm(KineticsPopulation,2,2); %p:2, distance to x axis or the origin in the high-dimensional space
% [MinIValue_vec,Min_idx] = min(D_vec);
% [MaxValue_vec,Max_idx ] = max(D_vec);
% para_MAX = ParaPopulation(  Max_idx,:  )'; %convert to vertical
% para_MIN = ParaPopulation(  Min_idx,:   )';
% C_MAX =  KineticsPopulation(Max_idx,:)' ; %vertical
% C_MIN =  KineticsPopulation(Min_idx ,:)';
end
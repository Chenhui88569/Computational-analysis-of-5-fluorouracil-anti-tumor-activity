function   [CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN] = Resistance_InferentialAnalysis(...
   KineticsPopulation, ParameterPopulation_size, ParameterPopulation     ) %CI: confidence interval

SEM_TS =std(KineticsPopulation )/sqrt(ParameterPopulation_size  ); %row vector
conf = 0.9;
alpha  = 1- conf;
pLo = alpha  /2;
pUp = 1- pLo ;
ts = tinv( [pLo  pUp ], ParameterPopulation_size-1 ); %90%confidence level 
c_mean = [ transpose(mean( KineticsPopulation))  transpose(mean( KineticsPopulation))];%the mean of each column
CI = c_mean + SEM_TS'*ts;  %time*2
C_mean = c_mean(:,1); %vertical
%find the overall max ;the maximum minimum
%zero crossing can not be applied to all the situations occurring in the
%model
MedianPopulation = median(  KineticsPopulation ,2  );%median of each row
[~ ,MedianIndex_inPopulation] = min( abs( MedianPopulation-median(MedianPopulation   )));%return the position of the
%first value that is closet to the median.
if MedianPopulation( MedianIndex_inPopulation   )>KineticsPopulation(  MedianIndex_inPopulation,1   ) %larger than initial value
    %upward curve , 
    MaxPoint = max(KineticsPopulation(:,2:end),[],2) ;%find the minimum value of each row;  except the initial time point
    Max_index  = find( MaxPoint == max(  MaxPoint  ) ) ; %maximum of the max
    Min_index = find( MaxPoint == min(  MaxPoint  ) );%maximum of the min 
elseif  MedianPopulation( MedianIndex_inPopulation   )< KineticsPopulation(  MedianIndex_inPopulation,1   )
     %downward curve ,
    MinPoint = min(KineticsPopulation(:,2:end),[],2) ;%find the minimum value of each row;  except the initial time point
    Max_index  = find( MinPoint ==max(  MinPoint  ) ) ; %maximum of the min
    Min_index = find( MinPoint == min(  MinPoint  ) );
end
%if the parameter sets corresponding to min are more than one
if size( Min_index,1) >1
     Min_index =  Min_index(1);
end
if size( Max_index,1) >1
     Max_index =  Max_index(1);
end
para_MAX = ParameterPopulation(  Max_index,:  )'; %convert to vertical
para_MIN = ParameterPopulation(  Min_index,:   )';
C_MAX =  KineticsPopulation(Max_index,:)' ; %vertical
C_MIN =  KineticsPopulation(Min_index ,:)';
end
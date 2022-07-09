function [idx,ClosetTimeCourse] = CompareTimeCourse(T_tumor_resistance,Cv_tumor_resistance ,...
    TumorVolume_PopulationKinetics,TumorVolume_TimeSpan)
    %  [dimension] note that two comparable time courses should be have the identical dimension
    %  That's why we need to use the interp1 to interpolate the elements of one of the time course at the quiry points.
    Cv_tumor_resistance_q =  interp1(T_tumor_resistance ,Cv_tumor_resistance,TumorVolume_TimeSpan' ,'PCHIP');
    Population_num = size(TumorVolume_PopulationKinetics,1);
    dist_col = zeros(Population_num,1);
     for i = 1: Population_num
        Kinetics_individual  = transpose( TumorVolume_PopulationKinetics(i,:) );
        %---------------------------------
        %  Hypthesis test of whether two set of time seris data are drawn from random samples from the normal distributions with the equal means.
        %  use [h,p] = ttest2(Cv_tumor_resistance_q , Kinetics_individual) can give rise to the same p value as the one way annov1 can.
        %  ttest2 is based on the two-sample t test 
        %  anova is based on the F-statistics
        %  h: the test decision for the null hypothesis with the possible values
        %            1: indicates the rejection of the null hypothesis at the alpha significance level, alpha  = 0.05
        %            0: indicates the fail to reject the null hypothesis 
        %--------------------------------------------
        %  F = anova1([ Cv_tumor_resistance_q Kinetics_individual  ] ,[],'off'); 
       % [h,p] = ttest2(Cv_tumor_resistance_q , Kinetics_individual )    ;
        dist = vecnorm(Cv_tumor_resistance_q -  Kinetics_individual,2  );
        dist_col(i)  =     dist;
     end
    [~,idx] = min(dist_col);
    ClosetTimeCourse =  TumorVolume_PopulationKinetics(idx,:);
  
end
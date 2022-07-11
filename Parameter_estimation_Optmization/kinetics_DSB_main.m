clear  
plasma_UnitConversion = 10^6/130.077; 
day2min  =60*24;
t_DSB_hour = [ 0
    8
    24
    32
    48
    64
    72];
%t_DSB_hour = t_DSB_hour- 24*ones(5,1);
t_DSB = t_DSB_hour*60;
t_DSB_xstick = reshape(t_DSB ,[1,size(t_DSB,1)]);
t_DSB_xlabel = reshape(t_DSB_hour ,[1,size(t_DSB,1)]) ;
t_DSB_xlabel  = num2cell(t_DSB_xlabel,1);


c_DSB = [42.755
    74.758
    74.276
    113.8229%121.733
    226.8139
    189.0788%163.286
    132.248
    ];



para_num_DSB = 7;
init_para = zeros(1 ,para_num_DSB);
%init_para =  reshape(init_para,[1,size(init_para,1)]);
para_totalnum = size(init_para,2);
population_size =100;
LB = 0.75;
UB = 1.25;


MaxIter =4;
Iter_num = 1;
theta_collection = zeros(MaxIter+1,para_num_DSB+1);
RSS_collection = zeros(MaxIter,1);
Optimizeoption_flag = 2;
multistart_Startpoint_flag =1; %use CustomizeStartPoint to genertate starting points1
myfun_DSB = @kinetics_DSB_withoutDDE;
myfun_DSB_plot = @kinetics_DSB_plot;

Dose_info = [
    45*day2min	3	7*day2min	 93/20*plasma_UnitConversion
    45*day2min	3	7*day2min	 93/10*plasma_UnitConversion
    45*day2min	3	7*day2min	 93/5*plasma_UnitConversion
    45*day2min	3	7*day2min	 93/2*plasma_UnitConversion
    45*day2min	3	7*day2min	 93*plasma_UnitConversion
    45*day2min	3	7*day2min	 0
    ];%Duration, Interval_TotalNum,DoseFrequency ,Dose
Dose_info (:,1:3)  = ones(6,1)* [t_DSB(end)+100*60   0   0	 ];
Dose_info_display = append("Model: ", ["5 mg/kg","10 mg/kg","20 mg/kg","50 mg/kg","100mg/kg"]);
 
estimation_flag =1; % change 
if estimation_flag==1
    if multistart_Startpoint_flag == 0 %RandomStart
        theta0 = rand(para_num_DSB,1);
        
        theta_collection(1,1:1:para_num_DSB) = theta0;
        disp('initial theta');
        disp(theta0);
        
        theta_collection(Iter_num,1:1:para_num_DSB ) = theta0;
        while Iter_num <= MaxIter
            start_point_flag = 0;%random start
            [theta_g,RSS] = kinetics5_para_estimation(t_DSB,c_DSB, theta0, myfun_DSB, para_num_DSB, Optimizeoption_flag, start_point_flag,[]) ;
            theta0 = theta_g;
            theta_collection(Iter_num+1,1:1:para_num_DSB) = theta_g;
            theta_collection(Iter_num+1,para_num_DSB+1) = RSS;
            RSS_collection(Iter_num) = RSS;
            fprintf('\t\tRSS in %dth iteration\t\t = %8.5f\n', Iter_num, RSS);
            Iter_num = Iter_num+1;
            
        end
        
    end
    
    if multistart_Startpoint_flag == 1 % CustomizeStart
        start_point_flag = 1;
        
        %theta0 =rand(1, para_num_DSB) ;
        %theta0 = ones(1, para_num_DSB);
        
        theta0 = [
%             30*60%20
%             0.71095
%             0.30610
%             0.11757
%             0.40340
2498.92633
2.01797
150.65770
10.46107
183.75192
0.03802
2.89273];
theta0 = theta0';

        init_para = theta0;
        theta_collection(Iter_num,1:1:para_num_DSB ) = theta0;
        c_DSB_weight = c_DSB;
        c_DSB_weight(c_DSB~=0) = 1;
        while Iter_num <= MaxIter
            %use the best fitted theta of the previous iteration as the initial value of parameters in the current iteration
            para_initGuess_modified = initGuess_generater(init_para,population_size,LB,UB,para_totalnum,[],[]);
            start_para = para_initGuess_modified;
            theta0 = init_para;
            [theta_g,RSS] = kinetics5_para_estimation(t_DSB, c_DSB_weight, theta0,myfun_DSB, para_num_DSB, Optimizeoption_flag, start_point_flag,start_para) ;
            init_para = theta_g; %1*para_num matrx
            theta_collection(Iter_num+1,1:1:para_num_DSB ) = theta_g;
            theta_collection(Iter_num+1,para_num_DSB +1) = RSS;
            RSS_collection(Iter_num) = RSS;
            for k1 = 1:length(theta_g)
                fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1,theta_g(k1))
            end
            fprintf('\t\tRSS in %dth iteration\t = %8.5f\n', Iter_num, RSS);
            pause(2);
            Iter_num = Iter_num+1;
            
             %plot
            f = figure(1);
            plot(t_DSB, c_DSB ,'o', 'MarkerSize',7,'DisplayName', 'data');%data
            hold on
            for i  = 1:1: size( Dose_info ,1 )-1
                [Tv,Cfit]  =  myfun_DSB_plot( theta_g ,  Dose_info(i,:));
                plot(Tv, Cfit, 'LineWidth',2, 'DisplayName',Dose_info_display(i)  );%model
            end
            hold off
            grid
            xticks(t_DSB_xstick)
            xticklabels(t_DSB_xlabel )
            xlabel('Time(h)')
            ylabel('γ-H2AX foci intensity(thousands)')
            title('Double strand break generation')
            legend
            drawnow;
            
            
            
        end
        
    end
    
  
    theta_fitted_DSB =  theta_collection ( find( RSS_collection == min( RSS_collection ))+1, 1:1:para_num_DSB);
    RSS_fitted_DSB = min(RSS_collection);
    for k1 = 1:length( theta_fitted_DSB)
        fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1,theta_fitted_DSB(k1))
    end
    fprintf('\t\tRSS\t\t = %8.5f\n',   RSS_fitted_DSB)
elseif estimation_flag ==0
    theta_fitted_DSB  = [
             2577.17019
10.50768
151.94553
10.06779
194.02545
0.12550
2.79398
% 1776.41779
% 0.69185
% 0.36925
% 0.11640
% 0.46572
        ];
    
end
%time curve

%plot
f = figure(1); 
plot(t_DSB(~ismember(t_DSB_hour, [32  64]) ), c_DSB(~ismember(t_DSB_hour, [32  64]))  ,'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'data');%data
hold on
plot([32 64]*60, c_DSB(ismember(t_DSB_hour, [32  64])) ,'*', 'MarkerSize',9,'DisplayName', 'added point');%data
hold on
for i  = 1:1: size( Dose_info ,1 )-1
    [Tv,Cfit]  =  myfun_DSB_plot(theta_fitted_DSB ,  Dose_info(i,:));
    plot(Tv, Cfit, 'LineWidth',2, 'DisplayName',Dose_info_display(i)  );%model
end
hold off
grid
xticks(t_DSB_xstick)
xticklabels(t_DSB_xlabel )
xlabel('Time(h)')
ylabel('γ-H2AX foci intensity(thousands)')
title('Double strand break generation')
legend


 










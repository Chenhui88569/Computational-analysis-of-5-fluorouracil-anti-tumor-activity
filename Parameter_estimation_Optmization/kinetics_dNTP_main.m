clear
Current_dNTP_M = 'alldNTP';
C_dNTP_measure_pertubation = dNTP_data_processing(Current_dNTP_M);
plasma_UnitConversion = 10^6/130.077; 
day2min = 24*60;
t_dNTP_hour =[
       0 
0.5 
1
2
4
6
8
9
12 
  ];
t_dNTP = t_dNTP_hour*60;
t_dNTP_xstick = reshape(t_dNTP ,[1,size(t_dNTP,1)]);
t_dNTP_xlabel = reshape(t_dNTP_hour ,[1,size(t_dNTP,1)]) ;
t_dNTP_xlabel  = num2cell(t_dNTP_xlabel,1);

 
 


para_num = 12;
init_para = zeros(1 ,para_num);
%init_para =  reshape(init_para,[1,size(init_para,1)]); 
population_size = 100;
LB = 0.75;
UB = 1.25;
Dose_info = [
    45*day2min	3	7*day2min	 93/20*plasma_UnitConversion
    45*day2min	3	7*day2min	 93/10*plasma_UnitConversion
    45*day2min	3	7*day2min	 93/5*plasma_UnitConversion
    45*day2min	3	7*day2min	 93/2*plasma_UnitConversion
    45*day2min	3	7*day2min	 93*plasma_UnitConversion
    45*day2min	3	7*day2min	 0
    ];%Duration, Interval_TotalNum,DoseFrequency ,Dose
Dose_info (:,1:3)  = ones(6,1)* [1*day2min  0   7*day2min	 ];
Dose_info_display = append("Model: ", ["5 mg/kg","10 mg/kg","20 mg/kg","50 mg/kg","100mg/kg"]);
addpoint = [2,4,8,9];


MaxIter =3;
Iter_num = 1;
theta_collection = zeros(MaxIter+1,para_num+1);
RSS_collection = zeros(MaxIter,1);
OptimizeSolver_flag = 2;%fmincon
multistart_Startpoint_flag =1; %use CustomizeStartPoint to genertate starting points1
myfun_dNTP = @kinetics_dNTP;
myfun_dNTP_plot = @kinetics_dNTP_plot;
 
estimation_flag = 1;
if estimation_flag ==1 
    if multistart_Startpoint_flag == 0 %RandomStart
        theta0 = rand(para_num,1);
        
        theta_collection(1,1:1:para_num) = theta0;
        disp('initial theta');
        disp(theta0);
        
        theta_collection(Iter_num,1:1:para_num ) = theta0;
        while Iter_num <= MaxIter
            start_point_flag = 0;%random start
            [theta_g,RSS] = kinetics5_para_estimation(t_DSB,c_DSB, theta0, myfun_dNTP, para_num, OptimizeSolver_flag, start_point_flag,[]) ;
            theta0 = theta_g;
            theta_collection(Iter_num+1,1:1:para_num) = theta_g;
            theta_collection(Iter_num+1,para_num+1) = RSS;
            RSS_collection(Iter_num) = RSS;
            fprintf('\t\tRSS in %dth iteration\t\t = %8.5f\n', Iter_num, RSS);
            Iter_num = Iter_num+1;
            
        end
        
    end
    
    if multistart_Startpoint_flag == 1 % CustomizeStart
        start_point_flag = 1;
       % theta0 =rand(1, para_num)*0.1 ;
         
  theta0 = [  
0.21351
32.50090
3.71407
73.34474
3.77027
0.30262
3.33454
1.35311
0.27614
35.66968
6.04530
0.07853
  ];
  theta0 = theta0';
init_para = theta0;
theta_collection(Iter_num,1:1:para_num ) = theta0;

    while Iter_num <= MaxIter
        %use the best fitted theta of the previous iteration as the initial value of parameters in the current iteration
        para_initGuess_modified = initGuess_generater(init_para,population_size,LB,UB,para_num,[],[]);
        start_para = para_initGuess_modified;
        [theta_g,RSS] = kinetics5_para_estimation(t_dNTP, C_dNTP_measure_pertubation, init_para,myfun_dNTP , para_num, OptimizeSolver_flag, start_point_flag,start_para) ;
        init_para = theta_g; %1*para_num matrx
        theta_collection(Iter_num+1,1:1:para_num ) = theta_g;
        RSS_collection(Iter_num) = RSS;
        for k1 = 1:length(theta_g)
            fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1,theta_g(k1))
        end
        fprintf('\t\tRSS in %dth iteration\t = %8.5f\n', Iter_num, RSS);
        pause(2)
        Iter_num = Iter_num+1;
        
        f = figure(1);
        data_plt=   plot(t_dNTP,  C_dNTP_measure_pertubation(:,1),'o', 'MarkerSize',8, 'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName',strcat('perturbation of ' , Current_dNTP_M,'-data' ) );%data
        hold on
        for i  = 1:1: size( Dose_info ,1 )-1
            [Tv,Cfit]  =  myfun_dNTP_plot( theta_g ,  Dose_info(i,:));
            plot(Tv, Cfit(:,1), 'LineWidth',2, 'DisplayName',Dose_info_display(i)  );%model  
        end
        hold off
        grid
        legend
        xticks(t_dNTP_xstick)
        xticklabels(t_dNTP_xlabel )
        xlabel('Time(h)')
        ylabel('pertubation')
        title(Current_dNTP_M)
        drawnow;
    end
    
    end
    
    
    theta_fitted_dNTP =  theta_collection ( find( RSS_collection == min( RSS_collection ))+1, 1:1:para_num);
    RSS_fitted_dNTP = min(RSS_collection);
    for k1 = 1:length( theta_fitted_dNTP)
        fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1,theta_fitted_dNTP(k1))
    end
    fprintf('\t\tRSS\t\t = %8.5f\n',  RSS_fitted_dNTP)
elseif estimation_flag == 0
    theta_fitted_dNTP = [
   0.21351
32.50090
3.71407
73.34474
3.77027
0.30262
3.33454
1.35311
0.27614
35.66968
6.04530
0.07853
% 1.36071
% 43.56677
% 28.79395
% 77.41199
% 4.12673
% 0.05359
% 2.08958
% 9.95190
% 0.85226
% 39.23792
% 6.54509
% 0.31116
];
end
 

f = figure(1);
C_data = C_dNTP_measure_pertubation(:,1);
plot( t_dNTP(~ismember(t_dNTP_hour,addpoint) ), C_data(~ismember(t_dNTP_hour,addpoint ) ),'o', 'MarkerSize',8,'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'data' );%data
hold on
plot(addpoint*60, C_data(ismember(t_dNTP_hour, addpoint ) ),'*', 'MarkerSize',9,'DisplayName', 'added point');%data
hold on
for i  = 1:1: size( Dose_info ,1 )-1
    [Tv,Cfit]  =  myfun_dNTP_plot(theta_fitted_dNTP,  Dose_info(i,:)); 
    hlp =  plot(Tv, Cfit(:,1), 'LineWidth',2, 'DisplayName',Dose_info_display(i)  );%model
 
end
hold off
grid
xticks(t_dNTP_xstick)
xticklabels(t_dNTP_xlabel )
xlabel('Time(h)')
ylabel('pertubation')
legend
title(Current_dNTP_M)








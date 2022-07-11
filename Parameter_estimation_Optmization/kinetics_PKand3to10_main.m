clear
model_flag  = 1;  % 1:pk; 2: compartment 3 to 9
Dataset_C = {};
Dataset_T = {};

%plasma
%t_plasma  = [0;7;15;30;40;50;60;70;80;90;105;120];
t_plasma  = [0;10;20;30;60;90;120;240];
c_plasma =  [ 92.9936306
70.2506231
67.0926138
48.329509
26.1156387 
11.3868339 
5.57423745
1.7525814]*10^6/130.077;
Dataset_C{1} = num2cell(c_plasma);
Dataset_T{1} = num2cell(t_plasma);
%  intracellular 5FU, FNUC
t_FNUC  = [0;15;30;40;50;60;70;80;90;105;120];
c_FNUC= [0     0 
159.9649225	29.33645133
133.3177433	50.21923414
110.8097048	42.81788951
39.85969015	49.19029524
86.67641041	43.56620871
62.37942122	75.70885706
11.87956738	67.8631979
46.23209588	63.11604794
6.243788366	73.77959661
8.477053493	61.79479684
];
Dataset_C{2} = num2cell(c_FNUC);
Dataset_T{2} = num2cell(t_FNUC);
  
%HMWA
c_RNA = [0
 9.61727039
6.487633256
9.155480343
2.365849466
2.880513238
0.338723071
0.041839875
] *10^3* 2*10^-3;
c_DNA = [0
   15.52683
13.606861
15.756227
9.8544088
14.428981
4.3345366
2.1127143  ]*0.2*10^-3; %
c_genome = [c_RNA c_DNA];
t_genome = [ 0; 2;8; 24;48;72;168;240  ];%h
t_genome = t_genome*60;%min 

Dataset_C{3} = num2cell( c_genome);
Dataset_T{3} = num2cell(t_genome); 

data_dUMP = [0  3.878
3 5.770
12 4.974
24 5.552
48 3.528] .*[60*ones(5,1)  ones(5,1)]; % time(h)  value

data_TS_complete  =[   0         0    0.0186
    2   0.0355         0
    6    0.0307         0
   24   0.0272    0.0135
   32    0.0217    0.0161
   48    0.0180    0.0228].*[60*ones(6,1)  ones(6,1) ones(6,1)]; %time(h) , TSFdUMP_complex , free TS




 

Dataset_C_g = Dataset_C;
Dataset_T_g =  Dataset_T;

if model_flag == 1 % plasma
    estimation_flag = 1; 
    T_data = cell2mat(Dataset_T{1});
    C_data =  cell2mat(Dataset_C{1});
   % [T_data_xtick, T_data_xlabel] = XtickLabelCreation( T_data , T_data );
    myfun = @kinetics_plasma;
    myfun_plot = @kinetics_plasma_plot;
    
    if  estimation_flag == 1
        MaxIter = 5;
        para_num  = 5;
        theta0 = rand(1,para_num);
       % norm_flag = 1;
        multistart_Startpoint_flag  = 1;  %use random StartPoint to genertate starting points1
        [theta_collection,RSS_collection,manytimes] = kinetics5_para_estimation( T_data, C_data, init_para,myfun, para_num, OptimizeSolver_flag, multistart_Startpoint_flag,start_para) ;

        [theta_collection,RSS_collection,manytimes ] = kinetics_iteration(...
            T_data, C_data, MaxIter, theta0 ,para_num , myfun, norm_flag , multistart_Startpoint_flag,0  );
        theta_fitted =  theta_collection ( find( RSS_collection == min( RSS_collection ))+1, 1:1:para_num);
        RSS_fitted = min(RSS_collection);
        for k1 = 1:length(theta_fitted)
            fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta_fitted(k1))
        end
        fprintf('RSS = %8.5f', RSS_fitted); 
        %plot
    elseif estimation_flag == 0
        theta_fitted = [10400.95432
5.66421
763.09044
625.02533
227261.94914];
    end
    tv = [min(T_data), max(T_data)+120*8];
    [Tv,Cfit]  =  myfun_plot( theta_fitted, tv);
    f = figure(1);
    Timespan = 0:120: max(T_data)+120*8;
    subplot(1,2,1)
    data_plt =  plot(T_data , C_data(:,1)/10^6*130.077 ,'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
    hold on
    plot(Tv, Cfit(:,1)/10^6*130.077 , 'LineWidth',2,'DisplayName','model');%model
    hold off 
    box off
    xlabel('Time(h)')
    ylabel('concentration(μg/mL)')
    xticks(Timespan)
    xticklabels(num2cell(Timespan/60))
    title('Central compartment')
    legend
    
    subplot(1,2,2)
%     data_plt =   plot( T_data , C_data(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
%     hold on
    plot(Tv, Cfit(:,2)/10^6*130.077, 'LineWidth',2,'DisplayName','model');%model
    hold off
    grid
    xlabel('Time(h)')
    ylabel('concentration(μg/mL)')
    xticks( Timespan)
    xticklabels(num2cell(Timespan/60))
    title('Peripheral compartment')
    legend
    %legend([hlp;data_plt], 'compartment1', 'compartment2', 'compartment1 data', 'compartment2 data','Location','SW')
    f = gcf;
    f.Position = [221.8,218.6,1134.4,512];
   
    f = figure(2);
    data_plt =  plot(T_data , C_data(:,1)/10^6*130.077 ,'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410]);%data
    hold on
    plot(Tv, Cfit/10^6*130.077 , 'LineWidth',2);%model
    hold off
    box off
    xlabel('Time(h)')
    ylabel('concentration(μg/mL)')
    xticks(Timespan)
    xticklabels(num2cell(Timespan/60))
    title('Central and peripheral compartment')
    legend({'data', 'central ','peripheral'})
    
end
population_size =200;
LB = 0.75; %0.75
UB = 1.25;%1.25
if model_flag == 2  % HWMA
    estimation_flag =0;  
    myfun = @kinetics_HMWA3to10_new;
    myfun_plot = @kinetics_HMWA3to10_plot;
    
    if  estimation_flag == 1
        MaxIter = 5;
        para_num  = 13; %21
        %theta0 =rand(1, para_num).*[5000*ones(1,2) ] ;
        theta0 = [
            %              1408.59448
            %             5415.53178
            %             rand(11,1)*0.1
            
%             1197.96765
%             5239.26919
%             rand(4,1)
%             rand(7,1)*0.01
1103.62127
4933.45078
0.01731
1.14684
0.00388
0.37143
4.67678
0.19896
0.03270
0.04440
0.10960
8.53048
6.39570            ];
        theta0 = theta0';
        
        norm_flag = 1;
        OptimizeSolver_flag = 2; %when the time frames of the species in a certain submodel are not consistent
        multistart_Startpoint_flag  = 1;  %use CustomizeStartPoint to genertate starting point
        Largerpertubation = 1;
        LargerParanum =2 ;
        init_para = theta0; %has to be row vector
        Iter_num = 1;
        theta_collection = zeros(MaxIter+1,para_num+1);
        RSS_collection = zeros(MaxIter,1);
        theta_collection(Iter_num,1:1:para_num ) = theta0;
        
        while Iter_num <= MaxIter
            %use the best fitted theta of the previous iteration as the initial value of parameters in the current iteration
            para_initGuess_modified = initGuess_generater(init_para,population_size,LB,UB,para_num, Largerpertubation, LargerParanum);
            start_para = para_initGuess_modified;
            [theta_g,RSS] = kinetics5_para_estimation([], [], init_para,myfun, para_num, OptimizeSolver_flag, multistart_Startpoint_flag,start_para) ;
            init_para = theta_g; %1*para_num matrx
             theta_collection(Iter_num+1,1:1:para_num ) = theta_g;
             theta_collection(Iter_num+1,para_num +1) = RSS;
             RSS_collection(Iter_num) = RSS;
             for k1 = 1:length(theta_g)
                 fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1,theta_g(k1))
             end
             fprintf('\t\tRSS in %dth iteration\t = %8.5f\n', Iter_num, RSS);
             pause(2);
             Iter_num = Iter_num+1;
             %plot
             [Tv_FNUC,Cfit_FNUC]  =  myfun_plot(theta_g, [min(t_FNUC), max(t_FNUC)]);
             f2=figure(2);
             subplot(1,3,1)
             hlp =  plot(Tv_FNUC, Cfit_FNUC(:,1), 'LineWidth',2, 'DisplayName','model');%model
             drawnow;
             grid
             xlabel('Time(min)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('5FU in interstitial fluid')
             legend('Location','northwest')
             
             subplot(1,3,2)
             data_plt =   plot(t_FNUC, c_FNUC(:,1),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             hold on
             hlp =  plot( Tv_FNUC, Cfit_FNUC(:,2), 'LineWidth',2, 'DisplayName','model');%model
             drawnow;
             hold off
             grid
             xlabel('Time(min)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('Intracellular 5FU')
             legend('Location','northwest')
             
             subplot(1,3,3)
             data_plt =   plot(t_FNUC, c_FNUC(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             hold on
             hlp =  plot(Tv_FNUC, Cfit_FNUC(:,3), 'LineWidth',2, 'DisplayName','model');%model
             drawnow;
             hold off
             grid
             xlabel('Time(min)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('5FU anabolites')
             legend('Location','northwest')
             f2= gcf;
             f2.Position = [-7,293,1539.2,420];
             
             
             [Tv_genome,Cfit_t_genome]  =  myfun_plot( theta_g, [min(t_genome), max(t_genome)]);
             
             figure(1)
             subplot(2,3,1)
             plot( t_genome,  c_genome(:,1),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             drawnow;
             hold on
             plot(Tv_genome, Cfit_t_genome(:,4), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xticks([24;48;72;168;240 ]*60)
             xticklabels({'1','2','3','7','10'} )
             xlabel('Time(day)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('F-RNA')
             legend
             
             subplot(2,3,2)
             data_plt=  plot(t_genome,  c_genome(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             drawnow;
             hold on
             hlp =  plot(Tv_genome, Cfit_t_genome(:,5), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xlabel('Time(day)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('F-DNA')
             xticks([24;48;72;168;240]*60)
             xticklabels({'1','2','3','7','10'} )
             legend
             
             subplot(2,3,3)
             data_plt=plot( data_TS_complete(:,1),  (data_TS_complete(:,2)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             drawnow;
             hold on
             hlp =  plot(Tv_genome, Cfit_t_genome(:,7), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xlabel('Time(h)')
             ylabel('the amount of activity(pmol/mg tissue)')
             xticks([24;48;72;168;240]*60)
             xticklabels({'1','2','3','7','10'} )
             title('TS-FdUMP')
             legend
             
%              subplot(2,3,4)
%              data_plt =   plot(data_TS_complete(:,1),  (data_TS_complete(:,3)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
%              drawnow;
%              hold on
%              hlp =  plot(Tv_genome, Cfit_t_genome(:,7), 'LineWidth',2, 'DisplayName','model');%model
%              hold off
%              grid
%              xlabel('Time(h)')
%              ylabel('the amount of activity(pmol/mg tissue)')
%              xticks([24;48;72;168;240]*60)
%              xticklabels({'1','2','3','7','10'} )
%              title('Free TS')
%              legend
             
             subplot(2,3,4)
             data_plt =   plot(data_dUMP(:,1),   data_dUMP(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','total dUMP data');%data
             drawnow;
             hold on
             hlp =  plot(Tv_genome, Cfit_t_genome(:, 6), 'LineWidth',2, 'DisplayName','totoal dUMP model');%model
             hold off
             grid
             xlabel('Time(h)')
             ylabel('the amount of activity(pmol/mg tissue)')
             xticks([24;48;72;168;240]*60)
             xticklabels({'1','2','3','7','10'} )
             title('Total dUMP')
             legend
             f= gcf;
             f.Position = [1,41,1536,746.4000000000001];
             
             [Tv_TS,Cfit_TS]  =  myfun_plot(theta_g,  min(t_genome):10:max( data_TS_complete(:,1)) );
             figure(3);
             subplot(1,3,1)
             data_plt=plot( data_TS_complete(:,1),  (data_TS_complete(:,2)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             drawnow;
             hold on
             hlp =  plot(Tv_TS, Cfit_TS(:,7), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xlabel('Time(h)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('TS-FdUMP')
             legend
%              
%              subplot(1,3,2)
%              data_plt =   plot(data_TS_complete(:,1),  (data_TS_complete(:,3)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
%              drawnow;
%              hold on
%              hlp =  plot(Tv_TS, Cfit_TS(:,7), 'LineWidth',2, 'DisplayName','model');%model
%              hold off
%              grid
%              xlabel('Time(h)')
%              ylabel('the amount of activity(pmol/mg tissue)')
%              title('Free TS')
%              legend
             
             subplot(1,3,2)
             plot(data_dUMP(:,1),   data_dUMP(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','total dUMP data');%data
             drawnow;
             hold on
             hlp =  plot(Tv_TS, Cfit_TS(:,6), 'LineWidth',2, 'DisplayName','totoal dUMP model');%model
             hold off
             grid
             xlabel('Time(h)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('Total dUMP')
             legend
             f3= gcf;
             f3.Position = [-7,293,1539.2,420];
             drawnow;
         end
         theta_fitted =  theta_collection ( find( RSS_collection == min( RSS_collection ))+1, 1:1:para_num);
         RSS_fitted = min(RSS_collection);
         for k1 = 1:length(theta_fitted)
             fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta_fitted(k1))
         end
         fprintf('RSS = %8.5f from %dth iteration ', RSS_fitted, find(RSS_collection==min(RSS_collection) ) );
    elseif estimation_flag == 0
        theta_fitted = [
        1064.24169	6537.154149	0.068161864	0.961440735	0.000612973	0.588098084	21.94138832	0.198449665	17.39699251	0.032339413	0.039127057	2.886080266	0.17527698	4.141404412

            ];
        %         1.41686
        % 29.60372
        % 1.14069
        % 11.97994
    end 
        
             [Tv_genome,Cfit_t_genome]  =  myfun_plot(theta_fitted, [min(t_genome), max(t_genome)]);
             
             figure 
             subplot(2,3,1)
             plot( t_genome,  c_genome(:,1),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data 
             hold on
             plot(Tv_genome, Cfit_t_genome(:,4), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xticks([24;48;72;168;240 ]*60)
             xticklabels({'1','2','3','7','10'} )
             xlabel('Time(day)')
             ylabel('amount (pmol/mg tissue)')
             title('F-RNA')
             legend
             
             subplot(2,3,2)
             data_plt=  plot(t_genome,  c_genome(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data 
             hold on
             hlp =  plot(Tv_genome, Cfit_t_genome(:,5), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xlabel('Time(day)')
             ylabel('amount (pmol/mg tissue)')
             title('F-DNA')
             xticks([24;48;72;168;240]*60)
             xticklabels({'1','2','3','7','10'} )
             legend
             
             subplot(2,3,3)
             data_plt=plot( data_TS_complete(:,1),  (data_TS_complete(:,2)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data 
             hold on
             hlp =  plot(Tv_genome, Cfit_t_genome(:,7), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xlabel('Time(h)')
             ylabel('amount (pmol/mg tissue)')
             xticks([24;48;72;168;240]*60)
             xticklabels({'1','2','3','7','10'} )
             title('TS-FdUMP')
             legend
             
             subplot(2,3,4)
             data_plt =   plot(data_TS_complete(:,1),  (data_TS_complete(:,3)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data 
             hold on
             hlp =  plot(Tv_genome, Cfit_t_genome(:,8), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xlabel('Time(h)')
             ylabel('amount (pmol/mg tissue)')
             xticks([24;48;72;168;240]*60)
             xticklabels({'1','2','3','7','10'} )
             title('Free TS')
             legend
             
             subplot(2,3,5)
              plot(data_dUMP(:,1),   data_dUMP(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','total dUMP data');%data 
             hold on
             hlp =  plot(Tv_genome, Cfit_t_genome(:,6), 'LineWidth',2, 'DisplayName','totoal dUMP model');%model
             hold off
             grid
             xlabel('Time(h)')
             ylabel('amount (pmol/mg tissue)')
             xticks([24;48;72;168;240]*60)
             xticklabels({'1','2','3','7','10'} )
             title('Total dUMP')
             legend
             f= gcf;
             f.Position = [1,41,1536,746.4000000000001];
             
              [Tv_FNUC,Cfit_FNUC]  =  myfun_plot(theta_fitted, [min(t_FNUC), max(t_FNUC)]); 
             figure;
             subplot(1,3,1)
             plot(Tv_FNUC, Cfit_FNUC(:,1), 'LineWidth',2, 'DisplayName','model');%model 
             grid
             xlabel('Time(min)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('5FU in interstitial fluid')
             legend('Location','northwest')
             
             subplot(1,3,2)
             plot(t_FNUC, c_FNUC(:,1),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             hold on
              plot( Tv_FNUC, Cfit_FNUC(:,2), 'LineWidth',2, 'DisplayName','model');%model 
             hold off
             grid
             xlabel('Time(min)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('Intracellular 5FU')
             legend('Location','northwest')
             
             subplot(1,3,3)
             plot(t_FNUC, c_FNUC(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             hold on
             plot(Tv_FNUC, Cfit_FNUC(:,3), 'LineWidth',2, 'DisplayName','model');%model 
             hold off
             grid
             xlabel('Time(min)')
             ylabel('the amount of activity(pmol/mg tissue)')
             title('5FU anabolites')
             legend('Location','northwest')
             f2 = gcf;
             f2.Position = [-7,293,1539.2,420]; 
             
             [Tv_TS,Cfit_TS]  =  myfun_plot(theta_fitted, [min(t_genome), max( data_TS_complete(:,1))]);
             [T_data_xtick_h, T_data_xlabel_h] = XtickLabelCreation( data_TS_complete(:,1),data_TS_complete(:,1)/60 );
             figure
             subplot(1,3,1)
             plot( data_TS_complete(:,1),  (data_TS_complete(:,2)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
  
             hold on
             plot(Tv_TS, Cfit_TS(:,7), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xticks(T_data_xtick_h)
             xticklabels(T_data_xlabel_h)
             xlabel('Time(h)')
             ylabel('amount (pmol/mg tissue)')
             title('TS-FdUMP')
             legend
             
             subplot(1,3,2)
             data_plt =   plot(data_TS_complete(:,1),  (data_TS_complete(:,3)),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             
             hold on
             hlp =  plot(Tv_TS, Cfit_TS(:,8), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xticks(T_data_xtick_h)
             xticklabels(T_data_xlabel_h)
             xlabel('Time(h)')
             ylabel('amount(pmol/mg tissue)')
             title('Free TS')
             legend
             
             subplot(1,3,3)
             plot(data_dUMP(:,1),   data_dUMP(:,2),'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
              
             hold on
             hlp =  plot(Tv_TS, Cfit_TS(:,6), 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xticks(T_data_xtick_h)
             xticklabels(T_data_xlabel_h)
             xlabel('Time(h)')
             ylabel('amount (pmol/mg tissue)')
             title('Total dUMP')
             legend
             f3= gcf;
             f3.Position = [-7,293,1539.2,420];
           
             
             figure(4)
             plot( data_TS_complete(:,1), data_TS_complete(:,3)*100 ./(data_TS_complete(:,3) + data_TS_complete(:,2)) ,'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName','data');%data
             
             hold on 
             plot(Tv_TS, Cfit_TS(:,8)*100./( Cfit_TS(:,8)+Cfit_TS(:,7) )  , 'LineWidth',2, 'DisplayName','model');%model
             hold off
             grid
             xticks(T_data_xtick_h)
             xticklabels(T_data_xlabel_h)
             xlabel('Time(h)')
             ylabel('Percentage')
             title(' %free TS')
             legend
end



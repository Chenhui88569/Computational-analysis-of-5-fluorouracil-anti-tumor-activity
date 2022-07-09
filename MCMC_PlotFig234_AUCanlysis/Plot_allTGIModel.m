%plot the simulation results of the nine models in one plot
%Each model occupy one subplot
%The information of all the models are returned by function @ CollectDataIntoNestedStructure
%The time course of tumor volume is determined by function  @Plot_kinetics_TumorVolume
% caculate the confidence interval of the time couses
function  Plot_allTGIModel(f2, theta_col, theta, Model_Index)
dbstop if error
newdir_parent  = 'Plot_TGI';
newdir = fullfile(newdir_parent, strcat("Model",num2str(Model_Index)) ) ;   % Your destination folder
label_add =  append ("(", ["A","B","C","D","E","F","G","H","I"] ,")");
AllModel_info   = CollectDataIntoNestedStructure ;
Color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.4940 0.1840 0.5560]];
Weektomin = 60*24*7;
day2min = 60*24;

TV_initial =  AllModel_info(Model_Index).InitialTv;
Dose_regime = AllModel_info(Model_Index).Dose_regime;
%Dose_regime_control =   Dose_regime ; Dose_regime_control(end) =0;
Dose = AllModel_info(Model_Index).FeaturedDose;
Dose_regime_cell = num2cell(Dose_regime);
ErrorBar = AllModel_info(Model_Index).ErrorBar ;
[Duration, Interval_TotalNum,DoseFrequency , ~ ] = Dose_regime_cell{:};
%label = ["Model "  num2str(Model_Index)  newline num2str(Dose) "mg/kg"  AllModel_info(Model_Index).Label];
label = [label_add(Model_Index) " Model "  num2str(Model_Index) ':' num2str(Dose) "mg/kg" newline  AllModel_info(Model_Index).Label];
label = strjoin(label); % to remove the line break.

Data_source = AllModel_info(Model_Index).DataSource;
Data = importdata( Data_source) ;
if Model_Index == 1
    Data(:,1) = Data(:,1)*day2min;
end

[Tsim,Ysim] = Plot_kinetics_TumorVolume(theta, Model_Index) ;
num_TimePoint = length(Tsim);
num_sample = 500;
Kinetics_Pop_Tv_treated = zeros(num_sample, num_TimePoint);
Kinetics_Pop_Tv_control = zeros(num_sample, num_TimePoint);
idx_shuffle = randperm(size(theta_col,1),num_sample);
theta_col = theta_col(idx_shuffle,:);

parfor m = 1: num_sample
    task  = getCurrentTask;
    taskid = task.Id;
    filename = ['logfile_' num2str(taskid) '.log'];
    outfile = fopen(  filename , 'a');
    theta_curr = theta_col(m,:);
    [Tfit_Tv ,Cfit_Tv] = Plot_kinetics_TumorVolume( theta_curr, Model_Index) ;
    %Cfit_col = Cfit_FNUC(:,[3,4,5])';
    
    Kinetics_Pop_Tv_treated(m,:) = Cfit_Tv(:,1)';  % Cfit_col(1,:) ;
    Kinetics_Pop_Tv_control(m,:) = Cfit_Tv(:,2)';  % Cfit_col(1,:) ;
    if rem(num_sample,50) == 0
        fprintf(  outfile ,'sample , %d \r', m);
    end
end

CI_treated = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_Tv_treated,  Ysim(:,1), num_sample    );

CI_control = TimeCourses_ConfidenceInterval(...
   Kinetics_Pop_Tv_control,  Ysim(:,2), num_sample    );


set(0, 'CurrentFigure', f2);
subplot(3,3,Model_Index)
p1 = plot(Tsim,Ysim(:,1) ,'-','LineWidth',1.5,'color', Color(1,:));  % treated group 
hold on 
[ha1 hb hc] = shadedplot(Tsim', transpose(CI_treated(:,1)) ,transpose(CI_treated(:,2)) , [17 17 17]/255,'none'); % confidence interval
hold on
p2 = plot(Tsim,Ysim(:,2) ,'-.', 'LineWidth',1.5, 'color',Color(2,:)) ; % control group;  different line style than the treated group
hold on
[ha2 hb hc] = shadedplot(Tsim', transpose(CI_control(:,1)) ,transpose(CI_control(:,2)) , [17 17 17]/255,'none'); % confidence interval
hold on
%set(f,  'position', [1,1,908,515])
%e1 = errorbar( Data(:,1) , Data(:,2),ErrorBar(:,1) , 'o','MarkerEdgeColor',Color(1,:));
e1 = errorbar( Data(:,1) , Data(:,2),ErrorBar(:,1),'o' );
%e1.Marker = 'o';
e1.MarkerSize = 10;
e1.Color = Color(1,:);
e1.CapSize = 15;
hold on
%errorbar( Data(:,1) , Data(:,3),ErrorBar(:,2), '^','MarkerEdgeColor',Color(2,:) ); %data
e2 = errorbar( Data(:,1) , Data(:,3),ErrorBar(:,2),'^'); %data
%e2.Marker = '^';
e2.MarkerSize = 10;
e2.Color = Color(2,:);
e2.CapSize = 15;
hold on
s = scatter( 0: DoseFrequency : Interval_TotalNum*DoseFrequency  , zeros( Interval_TotalNum+1,1   ), 80, 'filled','d' );
s.MarkerEdgeColor = Color(3,:);
s.MarkerFaceColor = Color(3,:);
hold off
box off
xticks( 0:Weektomin: fix( Duration /Weektomin)*Weektomin )
xticklabels(  num2cell( 0 :fix( Duration  /Weektomin )  ) )
xlabel('Time(weeks)')
ylabel(' Relative tumor volume')
title(label,'FontSize',14)
if Model_Index == 3
    legend([p1 p2 ha2 e1 e2 s],{'simulation-treated group', 'simulation-control group','',['Prediction interval' newline 'for control group'],...
    'data-treated group','data-control group','time of Injection'},...
    'Location','northeastoutside','FontSize', 14)
    legend boxoff 
end
% newdir_parent  = 'Plot_TGI';
% newdir = fullfile(newdir_parent, strcat("Model",num2str(Model_Index)) ) ;   % Your destination folder
% if ~isfolder( newdir )
%     mkdir( newdir  ) ;
% else
%     delete(  fullfile(newdir  ,'*' )   )
%     mkdir( newdir  ) ;
% end
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%     FigHandle = FigList(iFig);
%     FigName   = num2str(get(FigHandle, 'Number'));
%     set(0, 'CurrentFigure', FigHandle);
%     saveas( FigHandle, fullfile(newdir , [FigName '.png']));
% end
end
 function  f2 = PlotIntermediates_dNTP_DSB(f2, theta_col, theta)
Table_dir = 'DataSet_all.xlsx';
[T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir) ;
myfun_plot = @Plot_kinetics_dNTP_DSB;
[Tv_mean,Cfit_mean]  =  myfun_plot(theta);   %The difference in time span is already shown in the function inside;



num_sample = 500;
Kinetics_Pop_dNTP = zeros(num_sample, length(Cfit_mean{1}));
Kinetics_Pop_DSB = zeros(num_sample, length(Cfit_mean{2}));
idx_shuffle = randperm(size(theta_col,1),num_sample);
theta_col = theta_col(idx_shuffle,:);

dir_log = dir('*.log');
if ~isempty( [dir_log.name ])
    delete('*.log')
end

parfor i = 1:num_sample
    task  = getCurrentTask;
    taskid = task.Id;
    filename = ['logfile_' num2str(taskid) '.log'];
    outfile = fopen(  filename , 'a');
    theta_curr = theta_col(i,:);
    [~,Cfit_FNUC]  =  myfun_plot(theta_curr);
    %Cfit_col = Cfit_FNUC(:,[3,4,5])';
    Kinetics_Pop_dNTP(i,:) = Cfit_FNUC{1};  % Cfit_col(1,:) ;
    Kinetics_Pop_DSB(i,:) = Cfit_FNUC{2};
    if rem(num_sample,50) == 0
        fprintf(  outfile ,'sample , %d \r', i);
    end
end

label_add =  append ("(", ["A","B","C","D","E","F","G","I","H"] ,")");
set(0, 'CurrentFigure', f2);
set(f2, 'Position', get(0, 'Screensize'));
%% dNTP pool imbalance
subplot(3,3,8)
Anabolites_UnitConversion =10^6/130.077  ;
day2min = 60*24;
Dose_info = [7*day2min	 0	7*day2min	 92.9936306*Anabolites_UnitConversion]; % A week
t_dNTP =cell2mat( T_allsheets_output_cell(9)  ); 
t_dNTP_hour = t_dNTP/60; addpoint = [2,4,8,9];
[t_dNTP_xtick, t_dNTP_xlabel] = XtickLabelCreation(t_dNTP,t_dNTP_hour );
           
CI = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_dNTP,   Cfit_mean{1}', num_sample    ); % prediction interval
%f1 = figure;
Current_dNTP_M = 'alldNTP';
C_dNTP_measure_pertubation = dNTP_data_processing(Current_dNTP_M);
C_data = C_dNTP_measure_pertubation(:,1);
plot( t_dNTP(~ismember(t_dNTP_hour,addpoint) ), C_data(~ismember(t_dNTP_hour,addpoint ) ),'o', 'MarkerSize',6.5,'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'data' );%data
hold on
hlp =  plot(0:5:t_dNTP(end), Cfit_mean{1}, 'LineWidth',2 );%model
hold on
[ha hb hc] = shadedplot(0:5:t_dNTP(end), transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
hold off; box off 
xticks(t_dNTP_xtick)
xticklabels(t_dNTP_xlabel )
xlabel('Time(h)')
ylabel('pertubation')
title( strcat( label_add(8),  ' dNTP pool imbalance') ,'FontSize',15)
f1.CurrentAxes.TitleFontSizeMultiplier = 2;
newdir = 'Plot_Cellular_UpTodUMP';   % Your destination folder
%saveas( f, fullfile(newdir , ['dNTP' '.png']));

%%  DSB
t_DSB =cell2mat( T_allsheets_output_cell(10)  ); 
c_DSB = cell2mat( C_allsheets_output_cell(10)  ); 
t_DSB_hour = t_DSB/60 ;
[t_DSB_xtick, t_DSB_xlabel] = XtickLabelCreation(t_DSB,t_DSB_hour );

CI = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_DSB,   Cfit_mean{2}', num_sample    ); % prediction interval
%f2 = figure;
subplot(3,3,9)              
plot(t_DSB(~ismember(t_DSB_hour, [32  64]) ), c_DSB(~ismember(t_DSB_hour, [32  64]))  ,'o', 'MarkerSize',6.5,  'MarkerFaceColor',[0 0.4470 0.7410],'DisplayName', 'data');%data
hold on
plot(0:5:t_DSB(end), Cfit_mean{2}, 'LineWidth',2  );%model 
hold on
[ha hb hc] = shadedplot(0:5:t_DSB(end), transpose(CI(:,1)),transpose(CI(:,2)), [17 17 17]/255,'none'); % confidence interval
hold off; box off
xticks(t_DSB_xtick)
xticklabels(t_DSB_xlabel )
xlabel('Time(h)')
ylabel('γ-H2AX foci intensity(thousands)')
f2.CurrentAxes.TitleFontSizeMultiplier = 2;
title( strcat( label_add(9),  '  Double strand break generation'), 'FontSize',15)

% newdir = 'Plot_dNTP_DSB';   % Your destination folder
% if ~isfolder( newdir ) 
%    mkdir( newdir  ) ;
% else 
%   delete( fullfile(newdir,'*' )  ) ;
% end
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   %ohf = findobj( FigHandle); 
%   %ohf_axes = findobj(ohf , 'Type','axes' );
%   FigName   = num2str(get(FigHandle, 'Number'));
%   %set(0, 'CurrentFigure', FigHandle);
%   %set(FigHandle, 'Position', get(0, 'Screensize'));
%   saveas( FigHandle, fullfile(newdir , [FigName '.png']));
% end
end
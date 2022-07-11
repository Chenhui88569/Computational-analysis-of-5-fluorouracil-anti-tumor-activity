close all
clear
dbstop if error 
load Para_PK_col.mat Para_col S_curr
para_label = {'V_{max}' , 'Q_{21}' ,    'V_1',  'V_2',  'K_m','\sigma^2'};
Para_trun = Para_col(18000:end,:);
num_sample = size(Para_trun,1);
num_para = length(para_label);
Para_mean_col = zeros(1,5);
T_col =  cell(num_para ,1);

newdir = 'Plot_PK';   % Your destination folder 
if ~isfolder( newdir   ) 
   mkdir( newdir   ) ;
else 
  rmdir( newdir ,'s'  ) ;
  mkdir( newdir ) ;
end
f = figure;
set(f, 'Position', get(0, 'Screensize'));
for i = 1:6
    subplot(2,3,i)
    histogram(Para_trun(:,i) )
    hold on
    [counts, edges] = histcounts(Para_trun(:,i));
    locs = movmean(edges, 2, 'Endpoints', 'discard');
    plot(locs, counts, 'LineWidth', 3);
    ax = gca;
    title(para_label(i),'FontSize', 15)
    Para_mean = mean(Para_trun(:,i));
    Para_mean = round(Para_mean,3);
    Para_mean_col(i) =  Para_mean;
    hold on
    xline(Para_mean,'LineWidth',3)
    x = Para_trun(:,i);
    CI  = quantile(x, [0.25, 0.75]);
    %bar_label = strcat( '\bar', '{', para_label(i),'}'   );
    bar_label = strcat( '(',para_label(i),')','_{avg}' );
    temp = strcat( num2str( Para_mean) ,'(', num2str(CI(1)), ',', num2str(CI(2)),')');
    T_col{i} = temp;
    % text( ax.XTick(1),  (ax.YTick(end)+ax.YTick(end-1))/2  ,temp,'tex'  )
    hold off
end
sgtitle( 'Histograms for parameters in \Theta_1','FontSize', 15)
saveas( f, fullfile(newdir  , 'FigS1_S3.tif'));

f2 = figure;
clf
for i = 1:6
    subplot(2,3,i)
    plot(Para_col(:,i),'k-','LineWidth',1)
    xlabel('iteration')
    %axis([1+T*fix(curr_T/T),T+T*fix(curr_T/T),0,15]);
    grid off
    title(para_label(i),'FontSize',15)
end

c_plasma =  [ 92.9936306
70.2506231
67.0926138
48.329509
26.1156387 
11.3868339 
5.57423745
1.7525814]; % pmol/ML
t_plasma  = [0;10;20;30;60;90;120;240];
 
time_span_PK = 0:1:240;  num_TimePoint = length(time_span_PK);
[Tsim,Csim,C_diff] = kinetics_plasma(Para_mean_col,time_span_PK );
Kinetics_Pop_PK_central = zeros(num_sample, num_TimePoint);
Kinetics_Pop_PK_peripheral = zeros(num_sample, num_TimePoint);
parfor i = 1: num_sample
    task  = getCurrentTask;
    taskid = task.Id;
    filename = ['logfile_' num2str(taskid) '.log'];
    outfile = fopen(  filename , 'a');
    theta_curr = Para_trun(i,:);
    [~ ,Cfit,~] = kinetics_plasma( theta_curr ,time_span_PK) ;
    %Cfit_col = Cfit_FNUC(:,[3,4,5])';
    Kinetics_Pop_PK_central(i,:) = Cfit(:,1)';  % Cfit_col(1,:) ;
    Kinetics_Pop_PK_peripheral(i,:) = Cfit(:,2)';
    if rem(num_sample,100) == 0
        fprintf(  outfile ,'sample , %d \r', i);
    end
end

CI_central = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_PK_central, Csim(:,1), num_sample    ); % prediction interval
CI_peripheral = TimeCourses_ConfidenceInterval(...
  Kinetics_Pop_PK_peripheral, Csim(:,2), num_sample    ); % prediction interval
f_PK = figure;
d = plot(t_plasma  , c_plasma ,'o', 'MarkerSize',8,  'MarkerFaceColor',[0 0.4470 0.7410]);%data
hold on
f1 = plot(Tsim, Csim(:,1), 'LineWidth',2);%model; central and peripheral comparmtent are plotted together.
hold on
f2 = plot(Tsim, Csim(:,2), '-.','LineWidth',2);%model; central and peripheral comparmtent are plotted together.
hold on
f3 = plot([Tsim Tsim], CI_central, '--k','LineWidth',2);%model; central and peripheral comparmtent are plotted together.
hold on 
f4 = plot([Tsim Tsim], CI_peripheral, '.k','LineWidth',2);%model; central and peripheral comparmtent are plotted together.
% [ha1 hb hc] = shadedplot(Tsim', transpose(CI_central(:,1)),transpose(CI_central(:,2)), [17 17 17]/255,'none'); % confidence interval
% hold on
% [ha2 hb hc] = shadedplot(Tsim', transpose(CI_peripheral(:,1)),transpose(CI_peripheral(:,2)), [17 17 17]/255,'none'); % confidence interval
hold off
box off
xlabel('Time(min)')
ylabel('concentration(μg/mL)')
xticks(t_plasma)
xticklabels(num2cell(t_plasma))
title('Central and peripheral compartment','FontSize',15)
%legend([d f1 f2 ha1 ha2],{'Data:central', 'Simulation:central','simulation:peripheral','','Prediction level:central','','Prediction level:peripheral'},'FontSize',13)
legend([d f1 f2 f3(1) f4(1)],{'Data:central', 'Simulation:central','simulation:peripheral','Prediction level:central','Prediction level:peripheral'},'FontSize',13)
legend boxoff
%legend([d f1 f2],{'data', 'simulation-central','simulation-peripheral'},'location','northeastoutside','FontSize',15)
saveas( f_PK, fullfile(newdir, 'Fig2.tif'));
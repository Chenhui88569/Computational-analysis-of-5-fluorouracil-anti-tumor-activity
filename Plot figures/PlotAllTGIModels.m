%plot the simulation results of the nine models in one plot
%Each model occupy one subplot
%The information of all the models are returned by function @ CollectDataIntoNestedStructure
%The time course of tumor volume is determined by function @kinetics_TvPlot_DoseDifference
clear
close all
AllModel_info   = CollectDataIntoNestedStructure ;
Color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.4940 0.1840 0.5560]];
Weektomin = 60*24*7;
day2min = 60*24;
for i =1:9
    %subplot(3,3,i) %originally subplot
    figure
    theta  = AllModel_info(i).Para;
    TV_initial =  AllModel_info(i).InitialTv;
    Dose_regime = AllModel_info(i).Dose_regime;
    Dose_regime_control =   Dose_regime ; Dose_regime_control(end) =0; 
    Dose = AllModel_info(i).FeaturedDose;
    Dose_regime_cell = num2cell(Dose_regime);
    ErrorBar = AllModel_info(i).ErrorBar ;
    [Duration, Interval_TotalNum,DoseFrequency , ~ ] = Dose_regime_cell{:}; 
    label = ["Model "  num2str(i)  newline num2str(Dose) "mg/kg"  AllModel_info(i).Label]; 
    label = strjoin(label);
    Data_source = AllModel_info(i).DataSource;
    Data = importdata( Data_source) ;
    if i ==1
        Data(:,1) = Data(:,1)*day2min;
    end
    [Tsim,Ysim] = kinetics_TvPlot_DoseDifference(theta,TV_initial, Dose_regime) ;
    [Tsim_control,Ysim_control] = kinetics_TvPlot_DoseDifference(theta,TV_initial,Dose_regime_control) ; 
    plot(Tsim,Ysim(:,1)/TV_initial ,'LineWidth',2,'color', Color(1,:))
    hold on 
    plot(Tsim_control,Ysim_control(:,1)/TV_initial ,'LineWidth',2, 'color',Color(2,:))
    hold on 
    errorbar( Data(:,1) , Data(:,2),ErrorBar(:,1) , 'o','MarkerEdgeColor',Color(1,:))
    hold on
    errorbar( Data(:,1) , Data(:,3),ErrorBar(:,2), '^','MarkerEdgeColor',Color(2,:) ); %data
    hold on
    scatter( 0: DoseFrequency : Interval_TotalNum*DoseFrequency  , zeros( Interval_TotalNum+1,1   ),   '*'  ) 
    hold off
    box off
    xticks( 0:Weektomin: fix( Duration /Weektomin)*Weektomin )
    xticklabels(  num2cell( 0 :fix( Duration  /Weektomin )  ) )
    xlabel('Time(weeks)')
    ylabel(' Relative tumor volume')
    title(label,'FontSize',15)
end
%legend({'simulation-treated group', 'simulation-control group','data-treated group','data-control group','injection'})

newdir = fullfile('..','plot' , 'All9Model_plot' ) ;   % Your destination folder
newdir  =strcat(FolderName , ['/', 'q', num2str(idx)  ]   );
if ~isfolder( newdir ) 
   mkdir( newdir  ) ;
else
   delete(  fullfile(newdir  ,'*' )   )
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  saveas( FigHandle, fullfile(newdir , ['Model_', FigName '.png']));
end

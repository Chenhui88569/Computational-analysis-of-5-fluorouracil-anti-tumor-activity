function Dose_difference_AUCe_damage(AllModel_info )
% addpath("D:\pythonworkshop\matlab\PKPD_5FU\compartmental analysis\TimeCourses_DSB_Anabolites")
% AUCs are caculated over a week 
Anabolites_UnitConversion = 10^6/130.077; 
day2min = 24*60;
Regime_AWeek = [7*day2min	 0	7*day2min	 92.9936306*Anabolites_UnitConversion
    7*day2min	  3     day2min  93/2*Anabolites_UnitConversion
    7*day2min	  1*2-1	7*day2min/2 	 93/5*Anabolites_UnitConversion
    7*day2min	 1*3-1	7*day2min/3 	 93/5*Anabolites_UnitConversion
    7*day2min  7  day2min  93/20*Anabolites_UnitConversion
    7*day2min  1  7*day2min/2    93/5*Anabolites_UnitConversion
    7*day2min  3   48*60    93/5*Anabolites_UnitConversion
    7*day2min  3   2*day2min  93/2*Anabolites_UnitConversion
    7*day2min  3  2*day2min  93/2*Anabolites_UnitConversion
    
    ];
gamma_DSB   = 0.2;
gamma_hill = 0.2;
E_deviation =@(DSB_0, DSB) ( DSB -DSB_0 )./DSB_0  ; 

Cumulative_result = [];
Weektomin = 7*24*60;
Color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.4660 0.6740 0.1880]];

PlotFlag = 1;
f_a = figure;
f_DNA = figure;
for model_index = 1:size(Regime_AWeek,1) 
    %% Preparation
    Dose_regime = AllModel_info(model_index).Dose_regime;
    Para = AllModel_info(model_index).Para ;
    InitialTv = AllModel_info(model_index).InitialTv;
    Dose_regime_cell = num2cell(Dose_regime);
    [Duration, Interval_TotalNum,DoseFrequency , ~ ] = Dose_regime_cell{:};  
    para_individual = AllModel_info(model_index).Para;
    Dose =  AllModel_info(model_index).FeaturedDose;
    UpstreamKinetics = AllModel_info(model_index).UpstreamKinetics;
    load(UpstreamKinetics, 'T_sim_treated', 'Cv_anabolites_treated','Cv_DSB_treated','T_sim_control', 'Cv_anabolites_control','Cv_DSB_control' )
    AnabolitesDSB_TimeSpan_control = T_sim_control;
    DSB_Course_control = Cv_DSB_control;

    para_individual_cell  = num2cell( para_individual );
    [T_tumor, IC50,E_max_damage, EC50,lamda_g, p_max,lamda_d ]=   para_individual_cell{:};
    E_damage =@(DSB_deviation)   E_max_damage.*DSB_deviation.^gamma_DSB./(EC50.^gamma_DSB + DSB_deviation.^gamma_DSB    );
    E_anabolites = @(c_Anabolites)  IC50.^gamma_hill ./( IC50.^gamma_hill +c_Anabolites.^gamma_hill) ; %exposure-effect
    
    %% compute profile of intermediates within a week
    %to save time ;take advantage of the existing time courses
    AnabolitesDSB_TimeSpan  = T_sim_treated;
    Anabolites_Course  = Cv_anabolites_treated;
    DSB_Course  =  Cv_DSB_treated;
    Week_idx = find(     AnabolitesDSB_TimeSpan  <= 7*day2min) ;
    Week_idx = Week_idx(end);
    AnabolitesDSB_TimeSpan_AWeek  = AnabolitesDSB_TimeSpan(1:Week_idx );
    Anabolites_Course_AWeek =  Anabolites_Course(1:Week_idx );
    DSB_Course_AWeek =   DSB_Course(1:Week_idx );
    %DSB deviation
    DSB_Course_control_q_AWeek =  interp1(  AnabolitesDSB_TimeSpan_control  ,DSB_Course_control ,    AnabolitesDSB_TimeSpan_AWeek ,'PCHIP'  );
    DSB_deviation_AWeek = E_deviation(DSB_Course_control_q_AWeek,   DSB_Course_AWeek );
    E_damage_individual_AWeek = E_damage(DSB_deviation_AWeek); 
    E_Anabolites_individual_AWeek =   E_anabolites(Anabolites_Course_AWeek);
    
    %% profile of intermediates within a week over the entire course of treatment
    AnabolitesDSB_TimeSpan  = T_sim_treated;
    Anabolites_Course  = Cv_anabolites_treated;
    DSB_Course  =  Cv_DSB_treated;
    DSB_Course_control_q =  interp1(  AnabolitesDSB_TimeSpan_control  ,DSB_Course_control ,    AnabolitesDSB_TimeSpan ,'PCHIP'  );
    DSB_deviation  = E_deviation(DSB_Course_control_q,   DSB_Course  );
    E_damage_individual  = E_damage(DSB_deviation ); 
    E_Anabolites_individual  =   E_anabolites(Anabolites_Course );
   
    if PlotFlag == 1
        label = ["Model "  num2str(model_index)  newline num2str(Dose) "mg/kg"  AllModel_info(model_index).Label];
        label = strjoin(label);
        %% E_damage
        set(0, 'CurrentFigure', f_DNA)
        set(f_DNA, 'Position', get(0, 'Screensize'));
        subplot(3,3,model_index)
        [Tsim ,Ysim ] = kinetics_TvPlot_DoseDifference(para_individual,model_index) ;  %function
        Dose_regime_control = Dose_regime;
        Dose_regime_control(end) = 0;
        
        yyaxis left
        plot(Tsim ,Ysim(:,1),'LineWidth',2,'color',Color(1,:))
        hold on
        plot(Tsim,Ysim(:,2),'LineWidth',2,'color',Color(1,:))
        ylabel('Relative tumor volume change')
        yyaxis right
        plot(AnabolitesDSB_TimeSpan,  real(E_damage_individual),'LineWidth',2,'color',Color(2,:))
        hold on
        scatter( 0: DoseFrequency : Interval_TotalNum*DoseFrequency  , zeros( Interval_TotalNum+1,1   ),56,Color(3,:),'*' )
        hold off
        xticks( 0:Weektomin: fix( Duration /Weektomin)*Weektomin )
        xticklabels(  num2cell( 0 :fix( Duration  /Weektomin )  ) )
        xlabel('Time(weeks)')
        ylabel('Effect of DSB')
        title(label)
       % saveas(A, strcat(folderpath_DSB, "\Model",num2str(ModelIndex),".png"  ))
        
       %% E_anabolites
        set(0, 'CurrentFigure', f_a)
        set(f_a, 'Position', get(0, 'Screensize'));
        subplot(3,3,model_index)
        %B = figure('Position', [258.6,93,605.6,470.4]);
        yyaxis left
        plot(Tsim ,Ysim(:,1) ,'LineWidth',2,'color',Color(1,:))
        hold on
        plot(Tsim ,Ysim(:,2) ,'LineWidth',2,'color',Color(1,:))
        hold on
        scatter( 0: DoseFrequency : Interval_TotalNum*DoseFrequency  , zeros( Interval_TotalNum+1,1   ),56,Color(3,:),'*' )
        
        ylabel('Relative tumor volume change')
        yyaxis right
        plot(AnabolitesDSB_TimeSpan,  real(E_Anabolites_individual),'LineWidth',2,'color',Color(2,:))
        xticks( 0:Weektomin: fix( Duration /Weektomin)*Weektomin )
        xticklabels(  num2cell( 0 :fix( Duration  /Weektomin )  ) )
        xlabel('Time(weeks)')
        ylabel('Effect of Anabolites')
        title(label)
        %saveas(B, strcat(folderpath_Anabolites, "\Model",num2str(ModelIndex),".png"  ))
    end
    set(0, 'CurrentFigure', f_a)
    sgtitle('Tumor volume and E_{a}(t)')
    
    set(0, 'CurrentFigure', f_DNA)
    sgtitle('Tumor volume and E_{DNA}(t)')
%% calculate AUCe, DSB and AUCe, anabolites.
    AnabolitesDSB_TimeSpan_AWeek =  AnabolitesDSB_TimeSpan_AWeek./day2min;
    AUC_E_DSB = trapz( AnabolitesDSB_TimeSpan_AWeek , E_damage_individual_AWeek  );
    AUC_E_A_complement = trapz( AnabolitesDSB_TimeSpan_AWeek , E_Anabolites_individual_AWeek  );
    AUC_E_A =  AUC_E_A_complement ;%1*7 -AUC_E_Acomplement  
    AUC_DSB =  trapz( AnabolitesDSB_TimeSpan_AWeek ,  DSB_deviation_AWeek  );
    AUC_A =  trapz( AnabolitesDSB_TimeSpan_AWeek ,  Anabolites_Course_AWeek );
    Cumulative_result(:,model_index )= [real(AUC_E_A); real(AUC_E_DSB);real( AUC_A);real( AUC_DSB ) ];
    %Compute T>MIC
    % the reason of interpolating MIC is that the time point associated with MIC may
    % not be directly found out from the evaluation points. 
%     TimeAboveMIC_fun = @(T  ,C  ,idx, MIC)  interp1(C(idx:idx+1),T(idx:idx+1) , MIC   );  
%     MIC_DSB =  EC50;
%     find_MIC_idx_DSB  =  find( DSB_deviation >MIC_DSB );
%     if isempty(find_MIC_idx_DSB )
%         TimeAboveMIC_DSB  = 0;
%     else
%         MIC_idx_DSB = find_MIC_idx_DSB(end);
%         TimeAboveMIC_DSB = TimeAboveMIC_fun(AnabolitesDSB_TimeSpan , DSB_deviation,  MIC_idx_DSB ,MIC_DSB  );
%     end
%     
%     MIC_Anabolites = IC50;
%     find_MIC_idx_A  =  find(  Anabolites_Course >MIC_Anabolites);
%     if isempty(   find_MIC_idx_A)
%         TimeAboveMIC_Anabolites = 0;
%     else
%         MIC_idx_A = find_MIC_idx_A(end);
%         TimeAboveMIC_Anabolites = TimeAboveMIC_fun(AnabolitesDSB_TimeSpan , Anabolites_Course, MIC_idx_A ,MIC_Anabolites  );
%     end
%     disp(Label);
%    
%     fprintf ("\t\tAUCe of DSB: %8.8f ; AUCe of anabolites: %8.8f;AUC of DSB:  %8.8f; AUC of anabolites:  %8.8f ; \n ....",real(AUC_E_DSB),real(AUC_E_A),...
%        real( AUC_DSB ),real( AUC_A)  );
%    

end
%cellstr(string()) can convert a numeric array to a cell array of character
%vector.
%1*9 cell array of char is the same as 1*1 cell array
idx = 1:  size( Regime_AWeek,1 ) ;
idx = string(idx );
Variablelabel = append("model ",idx); 
Variablelabel  = cellstr(Variablelabel ); 
Cumulative_result_table  = array2table(Cumulative_result,  'VariableNames',Variablelabel ,...
    'RowNames',{'AUCe of anabolites','AUCe of DSB','AUC of anabolites','AUC of DSB'});
writetable(Cumulative_result_table, "AUC.xls")
disp(Cumulative_result_table  );
 

   % plot(AnabolitesDSB_TimeSpan  ,   E_Anabolites_individual ,  'DisplayName' ,strcat(Label ,"E_anabolites") )
 
newdir_1 = 'Effects_plot/DSB';   % Your destination folder
newdir_2 = 'Effects_plot/anabolites';
%newdir  =strcat(FolderName , ['/', 'q', num2str(idx)  ]   );
if ~isfolder( newdir_1 ) 
   mkdir( newdir_1  ) ;
else 
    delete( fullfile(   newdir_1 ,'*'   ));
end
if ~isfolder( newdir_2 ) 
   mkdir( newdir_2  ) ;
else 
    delete( fullfile(   newdir_2 ,'*'   ));
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigNum= get(FigHandle, 'Number');
  FigName   = num2str(FigNum );
  set(0, 'CurrentFigure', FigHandle);
  if  rem(FigNum,2 ) ==1
      FigName = ['DSB','_Model', num2str( ceil(  FigNum/2)   )];
      saveas( FigHandle, fullfile(newdir_1 , [FigName '.png']));
  elseif   rem(FigNum,2 ) == 0
      FigName = ['Ana','_Model', num2str(    ceil(  FigNum/2)    )];
      saveas( FigHandle, fullfile(newdir_2 , [FigName '.png']));
  end
end
   
end

 
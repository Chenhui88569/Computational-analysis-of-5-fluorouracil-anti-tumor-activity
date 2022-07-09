clear
close all
UnitConversion = [  10^6/130.077, 24*60,10^-3 ,7*24*60   ];
UnitConversion_cell = num2cell(UnitConversion  );
[plasma_UnitConversion, day2min  , mm3tocm3  , week2min  ] = UnitConversion_cell{:};
AllModel_info   = CollectDataIntoNestedStructure ;
% 100  mg/kg on Day1 of Weeks 1 to 6 of a eight-week cycle, for a total of 4 cycles
% 50 mg/kg twice a week of Weeks 1 to 6 of a eight-week cycle, for a total of 4 cycles
% 100  mg/kg on Day 1 of Weeks 1 to 3 followed by 50 mg/kg twice a week of Weeks 4 to 6, for a total of 4 cycles(32 weeks)
%  50 mg/kg twice a week of Weeks 1 to 3  followed by 100  mg/kg on Day 1 of Weeks 4 to 6, for a total of 4 cycles(32 weeks)
% Duration, number of Dose interval , dose frequency ,  Dose
% Same control group
Total_cases = 4;
f_T = figure;
set(f_T, 'Position', get(0, 'Screensize'));
Color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]]; %color 3: 100 mg/kg

for indi_case = 1:Total_cases
    set(0, 'CurrentFigure', f_T)
    subplot(2,2,indi_case)
    if indi_case  == 1
        Dose = 92.9936306; %92.9936306;
        DoseFrequency =  week2min ;    % for one cycle
        TreatmentDuration_PerCycle = 8*week2min;
        num_cycle = 4;
        Interval_TotalNum_PerCycle = 5;
        TotalDuration = TreatmentDuration_PerCycle * num_cycle ;
        Dose_info = [TotalDuration,  Interval_TotalNum_PerCycle ,DoseFrequency , TreatmentDuration_PerCycle,    num_cycle , Dose ] ;
        theta = AllModel_info(1).Para;
        TV_initial = AllModel_info(1).InitialTv;
        scatter_injection =  [];
        for i = 1: num_cycle
            scatter_injection = [scatter_injection  (i-1)*TreatmentDuration_PerCycle: DoseFrequency : (i-1)*TreatmentDuration_PerCycle + Interval_TotalNum_PerCycle*DoseFrequency ];
        end
        Dose_Order  = ones(1,length(scatter_injection));
        [Tsim,Ysim] = kinetics_TvPlot_DoseDifference(theta,Dose_info, TV_initial , scatter_injection, Dose_Order  );
        s = scatter(  scatter_injection , zeros( length( scatter_injection),1 ), 80, 'filled','d' );
        s.MarkerEdgeColor = Color(3,:);
        s.MarkerFaceColor = Color(3,:);
        hold on
        
        %% If we add the original dosing regimen
%         TreatmentDuration_PerCycle = 4*week2min;
%         num_cycle = 1;
%         Interval_TotalNum_PerCycle =  3;
%         Dose_info = [TotalDuration,  Interval_TotalNum_PerCycle ,DoseFrequency , TreatmentDuration_PerCycle,    num_cycle , Dose ] ;
%         [Tsim_org,Ysim_org] = kinetics_TvPlot_DoseDifference(theta,Dose_info, TV_initial  );
%         p3 = plot(Tsim_org,Ysim_org(:,1) ,'b--', 'LineWidth',1.5);  % treated group
        %legend({'Treated group under new regimen' , 'Control group','Time of Injection', 'Treated group under regimen 1'},'FontSize', 14)
        %legend boxoff
    elseif indi_case  == 2 % 50 mg/kg twice a week of Weeks 1 to 6 of a eight-week cycle, for a total of 4 cycles
        Dose = 93/2;
        DoseFrequency =  week2min/2 ;    % for one cycle
        TreatmentDuration_PerCycle = 8*week2min;
        num_cycle = 4;
        Interval_TotalNum_PerCycle = 6*2-1;
        TotalDuration = TreatmentDuration_PerCycle * num_cycle ;
        Dose_info = [TotalDuration,  Interval_TotalNum_PerCycle ,DoseFrequency , TreatmentDuration_PerCycle,    num_cycle , Dose ] ;
        theta_all = AllModel_info(2).Para;
        theta_control =  AllModel_info(1).Para;
        theta =  [theta_all(1:4)  theta_control(5:7)];
        TV_initial = AllModel_info(2).InitialTv;
        scatter_injection =  [];        

        for i = 1: num_cycle
            scatter_injection = [scatter_injection  (i-1)*TreatmentDuration_PerCycle: DoseFrequency : (i-1)*TreatmentDuration_PerCycle + Interval_TotalNum_PerCycle*DoseFrequency ];
        end
        Dose_Order  = ones(1,length(scatter_injection));
        [Tsim,Ysim] = kinetics_TvPlot_DoseDifference(theta,Dose_info, TV_initial  , scatter_injection,Dose_Order  );
        s = scatter(  scatter_injection , zeros( length( scatter_injection),1 ), 80, 'filled','o' );
        s.MarkerEdgeColor = Color(4,:); % 50 mg/kg
        s.MarkerFaceColor = Color(4,:);
        hold on
     
     elseif indi_case  == 3  %100  mg/kg on Day 1 of Weeks 1 to 3 followed by 50 mg/kg twice a week of Weeks 4 to 6, for a total of 4 cycles(24 weeks)
        %% first treatment in order
        Dose = 92.9936306; %92.9936306;
        DoseFrequency =  week2min ;    % for one cycle
        TreatmentDuration_PerCycle = 3*week2min;
        num_cycle = 4;
        Interval_TotalNum_PerCycle = 2;
        TotalDuration = TreatmentDuration_PerCycle * num_cycle ;
        theta_1 = AllModel_info(1).Para;
        TV_initial = AllModel_info(1).InitialTv;
        Dose_info_1 = [TotalDuration,  Interval_TotalNum_PerCycle ,DoseFrequency , TreatmentDuration_PerCycle,    num_cycle , Dose ] ;
        %% second treatment in order
        Dose = 93/2;
        DoseFrequency =  week2min/2 ;    % for one cycle
        TreatmentDuration_PerCycle = 5*week2min;
        num_cycle = 4;
        Interval_TotalNum_PerCycle = 3*2-1;
        TotalDuration = TreatmentDuration_PerCycle * num_cycle ;
        Dose_info_2 = [TotalDuration,  Interval_TotalNum_PerCycle ,DoseFrequency , TreatmentDuration_PerCycle,    num_cycle , Dose  ] ;
        theta_all = AllModel_info(2).Para;
        theta_control =  AllModel_info(1).Para;
        theta_2 =  [theta_all(1:4)  theta_control(5:7)];
        theta_col = [theta_1;
            theta_2];
        Dose_info = [  Dose_info_1 ;Dose_info_2 ];
        %calculation of injection time 
        TreatmentDuration_PerCycle_col  = [0; Dose_info(:,4)];
        DoseFrequency_PerCycle_col  =  Dose_info(:,3);
        Interval_TotalNum_PerCycle = Dose_info(:,2);
        num_cycle= Dose_info(1,5);
        TotalTreatmentDuration_PerCycle = sum(   TreatmentDuration_PerCycle_col );
        %% injection time
        scatter_injection = [];
        Treatment_dur_start = cumsum(TreatmentDuration_PerCycle_col);
        num_dose_PerCycle = size(Dose_info,1);
        Dose_Order = [];
        for indi_cycle = 1: num_cycle
            for indi_dose =1: num_dose_PerCycle
                indi_d = DoseFrequency_PerCycle_col(indi_dose );
                indi_interval =  Interval_TotalNum_PerCycle( indi_dose );
                scatter_injection = [ scatter_injection   (indi_cycle-1)* TotalTreatmentDuration_PerCycle  + Treatment_dur_start( indi_dose): indi_d  :...
                    (indi_cycle-1)* TotalTreatmentDuration_PerCycle +  Treatment_dur_start(indi_dose) +  indi_interval*indi_d  ];%+  indi_interval*indi_d ];
                Dose_Order = [Dose_Order    indi_dose*ones(1,indi_interval +1)  ];
            end
        end
        [Tsim,Ysim] = kinetics_TvPlot_DoseDifference( theta_col,Dose_info, TV_initial,scatter_injection, Dose_Order  );
        TotalDuration= sum(Dose_info(:,1));
        s1= scatter(  scatter_injection(    Dose_Order == 1) , zeros( length(scatter_injection(Dose_Order == 1)),1 ), 80, 'filled','d' );
        s1.MarkerEdgeColor = Color(3,:); % 100 mg/kg
        s1.MarkerFaceColor = Color(3,:);
        hold on
        s2= scatter(  scatter_injection(Dose_Order == 2) , zeros( length(scatter_injection(Dose_Order == 2)) ,1 ), 80, 'filled','o' );
        s2.MarkerEdgeColor = Color(4,:); % 50 mg/kg
        s2.MarkerFaceColor = Color(4,:);
 
    elseif indi_case  == 4
        %% first treatment in order
        Dose = 92.9936306; %92.9936306;
        DoseFrequency =  week2min ;    % for one cycle
        TreatmentDuration_PerCycle = 5*week2min;
        num_cycle = 4;
        Interval_TotalNum_PerCycle = 2;
        TotalDuration = TreatmentDuration_PerCycle * num_cycle ;
        theta_1 = AllModel_info(1).Para;
        TV_initial = AllModel_info(1).InitialTv;
        Dose_info_1 = [TotalDuration,  Interval_TotalNum_PerCycle ,DoseFrequency , TreatmentDuration_PerCycle,    num_cycle , Dose ] ;
        %% second treatment in order
        Dose = 93/2;
        DoseFrequency =  week2min/2 ;    % for one cycle
        TreatmentDuration_PerCycle = 3*week2min;
        num_cycle = 4;
        Interval_TotalNum_PerCycle = 3*2-1;
        TotalDuration = TreatmentDuration_PerCycle * num_cycle ;
        Dose_info_2 = [TotalDuration,  Interval_TotalNum_PerCycle ,DoseFrequency , TreatmentDuration_PerCycle,    num_cycle , Dose  ] ;
        theta_all = AllModel_info(2).Para;
        theta_control =  AllModel_info(1).Para;
        theta_2 =  [theta_all(1:4)  theta_control(5:7)];
        
        theta_col = [theta_2;
            theta_1];
        Dose_info = [  Dose_info_2 ;Dose_info_1 ];
        %calculation of injection time 
        TreatmentDuration_PerCycle_col  = [0; Dose_info(:,4)];
        DoseFrequency_PerCycle_col  =  Dose_info(:,3);
        Interval_TotalNum_PerCycle = Dose_info(:,2);
        num_cycle= Dose_info(1,5);
        TotalTreatmentDuration_PerCycle = sum(   TreatmentDuration_PerCycle_col );
        Dose_Order = [];
        %% injection time
        scatter_injection = [];
        Treatment_dur_start = cumsum(TreatmentDuration_PerCycle_col);
        num_dose_PerCycle = size(Dose_info,1);
        for indi_cycle = 1: num_cycle
            for indi_dose =1: num_dose_PerCycle
                indi_d = DoseFrequency_PerCycle_col(indi_dose );
                indi_interval =  Interval_TotalNum_PerCycle( indi_dose );
                scatter_injection = [ scatter_injection   (indi_cycle-1)* TotalTreatmentDuration_PerCycle  + Treatment_dur_start( indi_dose): indi_d  :...
                    (indi_cycle-1)* TotalTreatmentDuration_PerCycle +  Treatment_dur_start(indi_dose) +  indi_interval*indi_d  ];%+  indi_interval*indi_d ];
                Dose_Order = [Dose_Order    indi_dose*ones(1,indi_interval +1)  ];

            end
        end
        [Tsim,Ysim] = kinetics_TvPlot_DoseDifference( theta_col,Dose_info, TV_initial,scatter_injection,Dose_Order  );
        TotalDuration= sum(Dose_info(:,1));
        
        s2= scatter(   scatter_injection( Dose_Order == 2) , zeros( length(scatter_injection(Dose_Order == 2)),1 ), 80, 'filled','d' );
        s2.MarkerEdgeColor = Color(3,:); % 100 mg/kg
        s2.MarkerFaceColor = Color(3,:);
        hold on
        s1= scatter(   scatter_injection(Dose_Order == 1) , zeros( length(scatter_injection(Dose_Order == 1)) ,1 ), 80, 'filled','o' );
        s1.MarkerEdgeColor = Color(4,:); % 50 mg/kg
        s1.MarkerFaceColor = Color(4,:);    
    end
     p1 = plot(Tsim,Ysim(:,1) ,'-','LineWidth',1.5,'color', Color(1,:));  % treated group
     %[ha1 hb hc] = shadedplot(Tsim', transpose(CI_treated(:,1)) ,transpose(CI_treated(:,2)) , [17 17 17]/255,'none'); % confidence interval
     hold on
     %e_sim = errorbar(Tsim,Ysim(:,1),  abs(transpose(CI_treated (:,1)-Ysim(:,1)) ), abs(transpose(CI_treated(:,2)-Ysim(:,1))) ,'Color', Color(1,:),'Marker' ,'.')        ;
     p2 = plot(Tsim,Ysim(:,2) ,'-.', 'LineWidth',1.5, 'color',Color(2,:)) ; % control group;  different line style than the treated group
     box off
     xticks( 0:week2min: fix( TotalDuration /week2min)*week2min )
     xticklabels(  num2cell( 0 :fix( TotalDuration   /week2min )  ) )
     xlabel('Time(weeks)')
     ylabel(' Relative tumor volume')
     title(['Dosage regimen ' num2str(indi_case)] , 'FontSize', 15)
     if indi_case == 3
         legend({'Injection of 100 mg/kg ', 'Injection of 50 mg/kg' ,'Treated group-new regimen' , 'Control group'},'FontSize', 14,'Location','northeastoutside')
         legend boxoff
     end

end

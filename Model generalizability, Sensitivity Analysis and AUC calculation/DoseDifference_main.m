clear
close all
%addpath('D:\pythonworkshop\matlab\PKPD_5FU\compartmental analysis\DataSet_TGI')
addpath('DataSet_TGI')
UnitConversion = 10^6/130.077;
day2min = 24*60;
mm3tocm3 = 10^-3;
week2min = 7*24*60;
prompt = { '\fontsize{9} Do you want to compute AUCe or see the figures?(2:sensitivity analysis ;1: only compute AUC, 0:see the figures  )' };
dlgtitle = 'Dose difference flags input';
definput = {'1'};
dims = [2 30];
opts.Interpreter = 'tex';
opts.Resize = 'on';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts );
OnlyAUCe = str2double(answer);  
AllModel_info   = CollectDataIntoNestedStructure ;
switch OnlyAUCe
    case   2%sensitivity analysis 
    ComputeAUC = 2;
    FocusOnTumorVolume= 3; %empty
    DoseDifference_plot(AllModel_info,[],[],ComputeAUC,[],[], [],...
        [] ,FocusOnTumorVolume,  ....
        [],[],...
        [],[] ) ;
    case 1
    ComputeAUC = 1;
    FocusOnTumorVolume= 3;
    DoseDifference_plot(AllModel_info,[],[],ComputeAUC,[],[], [],...
        [] ,FocusOnTumorVolume,  ....
        [],[],...
        [],[] ) ;
    case 0
    prompt = {'\fontsize{9} Select the model parameters used in testing dose difference?(which seed dose is going to use: 5, 20mg/kg, 50 mg/kg or 100 mg/kg)',...
        ' \fontsize{9} Only concern about the tumor growth pattern and the schedules is applied accordingly(1: yes , 0: no)',...
       };
    dlgtitle = 'Dose difference flags input';
    definput = {'20','1'};
    dims = [2 30];
    opts.Interpreter = 'tex';
    opts.Resize = 'on';
    answer = inputdlg(prompt,dlgtitle,dims,definput,opts );
    ParaFromTreatedGroupData_Dose = str2double( answer{1}  );
    FocusOnTumorVolume =str2double( answer{2}  );
    %SameSchedule =str2double( answer{3}  );   %test
    
    switch ParaFromTreatedGroupData_Dose
        case 100
            Model_idx = 1;
            theta_Tv= AllModel_info(Model_idx).Para;
            TV_initial =  AllModel_info(Model_idx).InitialTv;
            Data_array = importdata(AllModel_info(Model_idx).DataSource );
            Data_array(:,1)  = Data_array(:,1)*day2min;
            
        case 50 %Twice a week
            Model_idx = 2;
            theta_Tv= AllModel_info(Model_idx).Para;
            TV_initial =  AllModel_info(Model_idx).InitialTv;
            Data_array = importdata(AllModel_info(Model_idx).DataSource );
        
        case 20
        prompt = ['\fontsize{9} Use which model?' newline '1 :The one fit to 20mg/kg,twice a week(HT29)' newline ...
            '2 :The one fit to 20mg/kg, three times a week(HT29)' newline ...
            '3: The one fit to 20mg/kg, twice a week(SW620) ' newline...
            '4: The one fit to 20mg/kg, every the other day(SW620) ) '];%specified as a character vector, cell array of character vectors, or string array
        dlgtitle = 'Dose difference flags input';
        definput = {'1'};
        dims = [2 30];
        opts.Interpreter = 'tex';
        opts.Resize = 'on';
        answer_model20 = inputdlg(prompt,dlgtitle,dims,definput,opts );
        answer_model20= str2double(answer_model20);
        switch answer_model20
            case  1
                Model_idx = 3;
                theta_Tv= AllModel_info(Model_idx).Para;
                TV_initial =  AllModel_info(Model_idx).InitialTv;
                Data_array = importdata(AllModel_info(Model_idx).DataSource ); 
                prompt = ['\fontsize{9} Do you want to use the \lambda_{g} and P_{max} of model 4 instead of model 3?(1 : yes; 0: no)' ];
                dlgtitle = 'Dose difference flags input';
                definput = {'1'};
                dims = [2 30];
                opts.Interpreter = 'tex';
                opts.Resize = 'on';
                answer_interior = inputdlg(prompt,dlgtitle,dims,definput,opts );
                answer_interior  = str2double(answer_interior{1});
                if answer_interior== 1 
                   Model_idx_ForControlConsistency = 4;
                   %lambda_g and pmax, lambda_d are cell line relevant ;
                   %swap them
                   ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                   theta_Tv(5:7) =  ControlGroupPara; 
                   TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv;
                end
            case 2
                Model_idx = 4;
                theta_Tv= AllModel_info(Model_idx).Para;
                TV_initial =  AllModel_info(Model_idx).InitialTv;
                Data_array = importdata(AllModel_info(Model_idx).DataSource ); 
                prompt = ['\fontsize{9} Do you want to use the \lambda_{g} and P_{max}  of model 3 instead of model 4?(1 : yes; 0: no)' ];
                dlgtitle = 'Dose difference flags input';
                definput = {'1'};
                dims = [2 30];
                opts.Interpreter = 'tex';
                opts.Resize = 'on';
                answer_interior = inputdlg(prompt,dlgtitle,dims,definput,opts ); 
                if str2double(answer_interior{1}) == 1
                   Model_idx_ForControlConsistency = 3; 
                   ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                   theta_Tv(5:7) =  ControlGroupPara;
                   TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv; 
               end
      
            case 3 %20, twice a week, SW620
                Model_idx = 6;
                theta_Tv= AllModel_info(Model_idx).Para;
                TV_initial =  AllModel_info(Model_idx).InitialTv;
                Data_array = importdata(AllModel_info(Model_idx).DataSource ); 
                prompt = ['\fontsize{9} Do you want to use the \lambda_{g} and P_{max}  of model 5 instead of model 6?(1 : yes; 0 : no)' newline ...
                    'Do you want to use the \lambda_{g} and P_{max}  of model 7 instead of model 6? (2 : yes; 0 : no) ' ];
                dlgtitle = 'Dose difference flags input';
                definput = {'1'};
                dims = [2 30];
                opts.Interpreter = 'tex';
                opts.Resize = 'on';
                answer_interior = inputdlg(prompt,dlgtitle,dims,definput,opts ); 
                if str2double(answer_interior{1}) == 1
                   Model_idx_ForControlConsistency = 5; 
                   ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                   theta_Tv(5:7) =  ControlGroupPara;
                   TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv;  
                elseif str2double(answer_interior{1}) == 2
                    Model_idx_ForControlConsistency = 7;
                    ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                    theta_Tv(5:7) =  ControlGroupPara;
                    TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv;
                end
            
            
            case 4 %20, every the other day 
                Model_idx =7;
                theta_Tv= AllModel_info(Model_idx).Para;
                TV_initial =  AllModel_info(Model_idx).InitialTv;
                Data_array = importdata(AllModel_info(Model_idx).DataSource );
                prompt = ['\fontsize{9} Do you want to use the \lambda_{g} and P_{max}  of model 7 instead of model 9?(1 : yes; 0 : no)' newline ...
                    'Do you want to use the \lambda_{g} and P_{max}  of model 8 instead of model 9? (2 : yes; 0 : no) ' ];
                dlgtitle = 'Dose difference flags input';
                definput = {'1'};
                dims = [2 30];
                opts.Interpreter = 'tex';
                opts.Resize = 'on';
                answer_interior = inputdlg(prompt,dlgtitle,dims,definput,opts );
            
            if str2double(answer_interior{1}) == 1
                    Model_idx_ForControlConsistency = 5;
                    ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                    theta_Tv(5:7) =  ControlGroupPara;
                    TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv;
            elseif str2double(answer_interior{1}) == 2
                    Model_idx_ForControlConsistency = 6;
                    ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                    theta_Tv(5:7) =  ControlGroupPara;
                    TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv;
            end
        end 
        case 5
            Model_idx = 5;
            theta_Tv= AllModel_info(Model_idx).Para;
            TV_initial =  AllModel_info(Model_idx).InitialTv;
            Data_array = importdata(AllModel_info(Model_idx).DataSource ); 
            prompt = ['\fontsize{9} Do you want to use the \lambda_{g} and P_{max}  of model 6instead of model 5?(1 : yes; 0: no)'  newline...
                '\fontsize{9} Do you want to use the \lambda_{g} and P_{max}  of model 7 instead of model 5?(2 : yes; 0: no)'  ];
            dlgtitle = 'Dose difference flags input';
            definput = {'1'};
            dims = [2 30];
            opts.Interpreter = 'tex';
            opts.Resize = 'on';
            answer_interior = inputdlg(prompt,dlgtitle,dims,definput,opts ); 
            if str2double(answer_interior{1}) == 1
                Model_idx_ForControlConsistency = 6;
                ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                theta_Tv(5:7) =  ControlGroupPara;
                TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv;
            elseif str2double(answer_interior{1}) == 2
                Model_idx_ForControlConsistency = 7;
                ControlGroupPara = AllModel_info(Model_idx_ForControlConsistency).Para(5:7);
                theta_Tv(5:7) =  ControlGroupPara;
                TV_initial = AllModel_info(Model_idx_ForControlConsistency).InitialTv;
            end
    end
    FeaturedDose = AllModel_info(Model_idx).FeaturedDose;
    theta_fitted_Tumor_ControlGroupPara =   theta_Tv ([end-2,end-1,end]);
    %Seed_schedule used in display
    Seed_Model_DisplayLabel =  AllModel_info(Model_idx).Label;
end

switch  FocusOnTumorVolume
    case 1
        %only consider the dynamics of tumor growth
        %The time evolution of tumor growth under different doses and the same
        %or different dose schedules are plot in the same figure if the flag
        %right below is set as 0. If the flag is set as 1, the tumor growth
        %under different doses and different schedules could be shown
        %in separate figures along with the corresponding data for
        %the purpose of validation.
        for i = 1:10
            if i >1
                %It doesn'y appear in the frist inquiry
                prompt = {'Do you want to continue? (1:yes, 0:no ))'};
                dlgtitle = 'inquiry';
                definput = {'1'};
                dims = [2 30];
                opts.Resize = 'on';
                answer = inputdlg(prompt,dlgtitle,dims,definput,opts );
                if str2double(answer) == 0
                    break
                end
            end
            switch ParaFromTreatedGroupData_Dose
                case 100
                    prompt = {'\fontsize{9} Do you want to see the dynamics of tumor volume in all the cases in separate figs(1: yes, 0 :no)',...
                        ' \fontsize{9} if yes, which dose regime is of your interest(50, 20, 10, 5)'};
                case 20
                    prompt = {'\fontsize{9} Do you want to see the dynamics of tumor volume in all the cases in separate figs(1: yes, 0 :no)',...
                        ' \fontsize{9} if yes, which dose regime is of your interest(100, 50, 20, 10, 5)'};
                case 50
                    prompt = {'\fontsize{9} Do you want to see the dynamics of tumor volume in all the cases in separate figs(1: yes, 0 :no)',...
                        ' \fontsize{9} if yes, which dose regime is of your interest(100, 20, 10, 5)'};
                case 5
                    prompt = {'\fontsize{9} Do you want to see the dynamics of tumor volume in all the cases in separate figs(1: yes, 0 :no)',...
                        ' \fontsize{9} if yes, which dose regime is of your interest(100, 20, 10, 5)'};
                    
            end
            dlgtitle = 'Dose difference flags input';
            definput = {'1','20'};
            dims = [2 30];
            opts.Resize = 'on';
            answer = inputdlg(prompt,dlgtitle,dims,definput,opts );
            
            FocusOnTumorVolume_SeparatePlots =  str2double( answer{1}  ); %0,1
            FocusOnTumorVolume_SeparatePlots_CurrentDose = str2double( answer{2}  ) ; % if FocusOnTumorVolume_SeparatePlots =1; 5, 10 ,20 ,100
            DoseDifference_plot([],TV_initial,Data_array,[],theta_fitted_Tumor_ControlGroupPara,FeaturedDose,   Seed_Model_DisplayLabel,...
                theta_Tv,FocusOnTumorVolume,  ....
                FocusOnTumorVolume_SeparatePlots,FocusOnTumorVolume_SeparatePlots_CurrentDose,...
                [],[] );
        end
        
    case 0
        prompt = {'\fontsize{9} Do you want to see all the intermediates in one figure(1: yes, 0 :no)'};
        dlgtitle = 'Dose difference flags input';
        definput = {'1'};
        dims = [2 30];
        opts.Resize = 'on';
        answer = inputdlg(prompt,dlgtitle,dims,definput,opts );
        IntermediateSubplot_InOneFig =  str2double( answer{1}  );
        if IntermediateSubplot_InOneFig == 0
            %Plot the time profiles of the essential intermediates in separate figures
            %In each figure, the time profiles of the sepecies under different
            %dosage are plotted. The duration of model simulation can be
            %close to the one of experimental measurements or be fixed as two
            %weeks.
            TruncateDuration_ToDataTimeFrame = 1;
            DoseDifference_plot([],TV_initial,Data_array,[],theta_fitted_Tumor_ControlGroupPara,FeaturedDose,  Seed_Model_DisplayLabel,...
                theta_Tv,FocusOnTumorVolume,  ....
                [],[],...
                IntermediateSubplot_InOneFig,TruncateDuration_ToDataTimeFrame ) ;
        else
            %Plot the time profiles of the essential intermediates in one figure.
            %every time profile of the intermediate occupies one subplot.
            DoseDifference_plot([],TV_initial,Data_array,[],theta_fitted_Tumor_ControlGroupPara,FeaturedDose,  Seed_Model_DisplayLabel,...
                theta_Tv,FocusOnTumorVolume,  ....
                [],[],...
                IntermediateSubplot_InOneFig,[] ) ;
        end
end


function DoseDifference_plot(AllModel_info,TV_initial,Data_array,ComputeAUCe,theta_fitted_Tumor_ControlGroupPara,FeaturedDose,  Seed_Model_DisplayLabel,...
    theta_Tv,FocusOnTumorVolume,  ....
   FocusOnTumorVolume_SeparatePlots,FocusOnTumorVolume_SeparatePlots_CurrentDose,...
   IntermediateSubplot_InOneFig,TruncateDuration_ToDataTimeFrame    )

if ComputeAUCe == 1
    Dose_difference_AUCe_damage(AllModel_info);
%     [Dose_info_uniq_AUC,~,~] = unique(Dose_info,'rows','stable');
%     AbsoluteFrequency =fix(7*day2min./ Dose_info_uniq_AUC(1:end-1,3)) ;
%     DoseInfo_AUC = Dose_info_uniq_AUC(1:end-1,:); %no control
%     
%         AUC_Interstitial_col = [];
%         AUC_Edrug_col = [];
%         DurationInweek = 2;
%     DoseInfo_AUC(:,1) = DurationInweek*week2min; %Duration are the same
%     DoseInfo_AUC(:,2) = DurationInweek.*AbsoluteFrequency -1;
%     DoseInfo_AUC(4,2) =  DoseInfo_AUC(4,2)+1;
%     E_drug = @(c_Interstitial) E_max.*c_Interstitial.^3./(  EC50.^3 +c_Interstitial.^3   );
%     fprintf( strcat("Model parameters based on ", num2str( FeaturedDose ), " mg/kg dosage are used\n" ) );
%     for i =  1: size(DoseInfo_AUC,1)  %5 mg/kg, 10 mg/kg, 20mg/kg,50,100mg/kg,
%         fprintf( " \t The dosage model simulates: %s\n", Dose_info_display(i))  ;
%         [T_up2DSB, C_up2DSB] = kinetics_uptoDSB( DoseInfo_AUC(i,:) );
%         C_up2DSB_Interstitial  = C_up2DSB(:,1)/UnitConversion; %mug /mL
%         T_up2DSB_day = T_up2DSB/day2min;
%         AUC_Interstitial = trapz( T_up2DSB_day  , C_up2DSB_Interstitial  ); %Trapezoidal numerical integration
%         AUC_Interstitial_col (i) = AUC_Interstitial;
%         fprintf(" \t\t Area under the blood (or Interstitial) concentration–time curve over a period of  " +DurationInweek + "week(s)"  +"  = %8.5f\n", AUC_Interstitial  );
%         
%         AUC_Edrug = trapz(T_up2DSB_day ,  E_drug( C_up2DSB_Interstitial ) );
%         AUC_Edrug_col(i) = AUC_Edrug;
%         fprintf(" \t\t Area under the E_max curve over a period of  " +DurationInweek + "week(s)" + "  = %8.5f\n",AUC_Edrug);
%     end
%     disp(AUC_Interstitial_col')
%     disp(AUC_Edrug_col')
else
    
    t_TV_array = Data_array(:,1) ;
    Relative_data_ControlGroup = Data_array(:,3);
    
    theta_ControlGroup = theta_fitted_Tumor_ControlGroupPara;
    UnitConversion = 10^6/130.077;
    day2min = 24*60;
    mm3tocm3 = 10^-3;
    week2min = 7*24*60;
    TableProp_NormalizeTV_flag =1; %1:divide the tumor volume measurements by the starting value
    parameter_label_Control = ["\lambda_{g}","p_{max}","\lambda_{d}"];
    parameter_label = ["T_{tumor}","IC_{50}", "E_{max,damage}","EC_{50,damage}", parameter_label_Control ];
    parameter_unit_label = [" min", " pmol/mg" ," units"," units", "min^{-1}", " cm^3", " min^{-1}"]; %total 6 parameters
    
    %RNA_UnitConversion = 10^3* 2*10^-3;
    %DNA_UnitConversion = 0.2*10^-3;
    Table_dir = 'DataSet_all.xlsx';
    [T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir) ;
    
    
    %%%%
    myfun_kinetics  = @kinetics_TvPlot_DoseDifference; %kinetics_TvPlot_DoseDifference(theta,TV_initial, Dose_info)
    %Dose_label no 50 mh/kg
    DoseLabel = ["5 mg/kg"; "10 mg/kg" ; "20 mg/kg";"100 mg/kg"];%string array, each element is specified as double quote
    dose_c0 = [ 1/20; 1/10;1/5;1 ]* 92.9936306 *UnitConversion;
    
    t_Interstitial_5and10and20 = [0	0	 0
        0.08	0.08	0.08
        0.15	0.25	0.25
        0.25	0.5	0.5
        0.5	1	1
        0.75	1.5	1.5
        1	2	2
        ] *60; 
    c_Interstitial_5and10and20  = [ 14.55  30.82  60.38;
        4.51	12.07	31.39
        2.15	3.63	9.71
        1.38	1.59	4.38
        0.46	0.38	1.36
        0.15	0.13	0.58
        0.02	0.06	0.26
        ] ;
    
    %The reason of keep the original Dose_info matirx the same as the
    %varient one is that the command lines below shoould be common for both
    %of the cases to maintain the simplicity of the codes.
    Dose_info = [
        %  5*7*day2min 	5*3-1	7*day2min/3	 92.9936306*UnitConversion/20
        9*day2min 	8	 day2min 	 92.9936306*UnitConversion/20
        3*7*day2min	3*3-1	 7*day2min/3	 30.82*UnitConversion
        5*7*day2min	5*2-1	7*day2min/2	 30.82*UnitConversion
        2*7*day2min    6     2*day2min    60.38*UnitConversion
        3*7*day2min	 3*2-1	7*day2min/2	 60.38*UnitConversion
        7*day2min       3          day2min           92.9936306*UnitConversion/2
        3*7*day2min	 3*3-1	7*day2min/3	 60.38*UnitConversion
        45*day2min	     3	        7*day2min	      92.9936306*UnitConversion
        45*day2min	     3	        7*day2min	     0
        ];%Duration, Interval_TotalNum,DoseFrequency ,Dose4
    Dose_info_display = append( "Model:", [
        "5 mg/kg, daily for 9 days" ,...
        "10 mg/kg, 3 times a week for 3 weeks" , "10 mg/kg, twice a week for 5 weeks", ...
        "20 mg/kg, every other day for 2 weeks" , "20 mg/kg, twice a week for 3 weeks",...
        "50 mg/kg, daily, 4 times",...
        "20 mg/kg, 3 times a week for 3 weeks",...
        "100mg/kg, Weekly schedule for 4 weeks" ," control group" ]) ;
    Dose_info(:,end) = [ 1/20; 1/10;1/10;1/5 ;1/5; 1/2;1/5;1 ;0 ]*Dose_info(end-1,end);
    
    line_style = ["-", "-", "--", "--", ":", ":","-.","-." ];
    MarkerEdgeColor = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.4940 0.1840 0.5560]];
    Elongation = 0*week2min; %4
    
end

if FocusOnTumorVolume ==1  %focus on tumor volume
    if FocusOnTumorVolume_SeparatePlots == 0
        [Dose_info_uniq,~,~] = unique(Dose_info,'rows','stable');
        %If the schedules are the same when it come to the different
        %doses, the repetitions of the dose information would influence
        %the allocation of legend to plotted time series due to the
        %inconsistency in the dimension of ' Dose_info_display' string
        %vector and 'Dose_info' matrix. Therefore, the repetitive rows
        %should be removed, and at the same time the oreder specified in the
        %original matrix should be retained.
        Dose_info_uniq(:,1) = (Dose_info_uniq(end,1)+ Elongation )*ones(size(  Dose_info_uniq(:,1)  ) ); %Duration are the same
        [Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup] =  myfun_kinetics(theta_Tv, TV_initial ,Dose_info_uniq(end,:) );
        figure
        plot( t_TV_array, Relative_data_ControlGroup,'^' ,'MarkerSize',8    ) %control group
        hold on
        for i =  1: size(Dose_info_uniq,1)-1  %5 mg/kg, 10 mg/kg, 20mg/kg,100mg/kg,
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_uniq(i,:) );
            plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial ,line_style(i), 'LineWidth',1.5 );%model
        end
        hold on
        plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial ,line_style(end), 'LineWidth',1.5 );%ControlGroup use separate set of parameters
        hold off
        grid
        xlabel('Time(weeks)')
        xticks( 0:week2min: fix(   Dose_info_uniq(end,1) /week2min)*week2min )
        xticklabels(  num2cell( 0 :fix( Dose_info_uniq(end,1)/week2min )  ) )
        ylabel('normalized mean tumor volume (fold change)')
        legend([ "data control group"+" from "+num2str(FeaturedDose)+"mg/kg" + Seed_Model_DisplayLabel...
            Dose_info_display] ,'Location', 'best','FontSize',12) %string array
        title("Model simulation of Tumor Voume" + newline+...
            "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
            "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
            "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup) )), 'FontSize',12)
        
    elseif FocusOnTumorVolume_SeparatePlots == 1
        Dose_info_elongated = Dose_info;
       % Dose_info_elongated(:,1) = (Dose_info(end,1)+ Elongation )*ones(size(  Dose_info(:,1)  ) ); %no uniq
       Dose_info_elongated(:,1) = (2*week2min+ Elongation )*ones(size(  Dose_info(:,1)  ) ); %no uniq
        [Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup] =  myfun_kinetics(theta_Tv, TV_initial, Dose_info_elongated(end,:) )  ;
        if FocusOnTumorVolume_SeparatePlots_CurrentDose == 20
            CellLine =  ["SW620","HT29","SW480\beta 6"];
            Data_20mgkg_TwiceAWeek_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'20mgkg_TwiceAweek.xls' );
            Sheets_Dataset_Num = size(Data_20mgkg_TwiceAWeek_cell,1);
            %Individual growth pattern of the control group are fixed first 
            theta_Tv_data1 = [ 0.000230002 13.204156452 0.000005769]; %SW620
            theta_Tv_data2 = [0.00023 1.3684  0.00006];   %HT29 
            theta_Tv_data3  = [0.23844/day2min  0.99672   6.5382e-05/day2min]; 
            TV_initial_col = [ 0.0565    100*mm3tocm3  100*mm3tocm3];
            [Tv_Tumor_data2,Cfit_Tumor_variation]  = myfun_kinetics([theta_Tv(1:4) theta_Tv_data2] ,TV_initial_col(2), Dose_info_elongated(5,:)); %20 mg/kg
            [Tv_Tumor_data1,Cfit_Tumor_data1]  = myfun_kinetics([theta_Tv(1:4) theta_Tv_data1]     ,TV_initial_col(1), Dose_info_elongated(5,:)); %20 mg/kg
            [Tv_Tumor_data3,Cfit_Tumor_data3]  = myfun_kinetics([theta_Tv(1:4) theta_Tv_data3] ,TV_initial_col(3), Dose_info_elongated(5,:)); %20 mg/kg
            %data1 from SW620,we need to generalize the model on this cell line 
            figure('Position', [258.6,93,605.6,470.4]); 
            %Data points are normalized to the their own starting value
            %Simulated curves are normalized to their own initial
            %condition 
  
            if contains( Seed_Model_DisplayLabel, 'SW620')  %main 5 to 6
                 % if the dynamics of the control group of model 8 need to be reproduced.
                 % mix the dose related part of the model that is need to be validated and the natural growth part of the target model
                %The elements that the `theta_ControlGroup` contains could be either the control group parameters of the source model 
                %or the control group of the target model that need to be  captured
               plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial ,'--' ,'LineWidth',2,'color', [0.4660 0.6740 0.1880] );%model in the absence of treatment
               hold on
               [Tv_Tumor_data1,Cfit_Tumor_data1]  = myfun_kinetics( theta_Tv ,TV_initial, Dose_info_elongated(5,:)); %20 mg/kg
               Data_20mgkg_TwiceAWeek = cell2mat( Data_20mgkg_TwiceAWeek_cell(1));
               errorbar( Data_20mgkg_TwiceAWeek(:,1), Data_20mgkg_TwiceAWeek(:,2) ,Data_20mgkg_TwiceAWeek(:,4),'^','MarkerSize',8  ,'MarkerEdgeColor',MarkerEdgeColor(1,:) );
               hold on
               errorbar( Data_20mgkg_TwiceAWeek(:,1), Data_20mgkg_TwiceAWeek(:,3) ,Data_20mgkg_TwiceAWeek(:,5) ,'*','MarkerSize',8,'MarkerEdgeColor',MarkerEdgeColor(1,:) );
               hold on 
               plot(Tv_Tumor_data1,  Cfit_Tumor_data1(:,1)/TV_initial , 'LineWidth',2, 'color',MarkerEdgeColor(1,:) );%model
               hold off
                legend( ["model: no treatment" , ...
                strcat("data-", " Protocol 6", ",Cell line ",CellLine(1) ),...
                strcat("data-", " control group") ,...
                "Model simulation"  
                ], ...
                'Location', 'best', 'FontSize',12) 

            %data1 from HT29 ,we need to generalize the model on this cell line 
            elseif contains( Seed_Model_DisplayLabel, 'HT29')
               [Tv_Tumor_data2,Cfit_Tumor_data2]  = myfun_kinetics( theta_Tv ,TV_initial, Dose_info_elongated(5,:)); %20 mg/kg
               Data_20mgkg_TwiceAWeek = cell2mat( Data_20mgkg_TwiceAWeek_cell(2));
               errorbar( Data_20mgkg_TwiceAWeek(:,1), Data_20mgkg_TwiceAWeek(:,2) ,Data_20mgkg_TwiceAWeek(:,4),'^','MarkerSize',8  ,'MarkerEdgeColor',MarkerEdgeColor(2,:) );
               hold on
               errorbar( Data_20mgkg_TwiceAWeek(:,1), Data_20mgkg_TwiceAWeek(:,3) ,Data_20mgkg_TwiceAWeek(:,5),'*','MarkerSize',8,'MarkerEdgeColor',MarkerEdgeColor(2,:) );
               hold on 
               plot(Tv_Tumor_data2,  Cfit_Tumor_data2(:,1)/TV_initial , 'LineWidth',2, 'color',MarkerEdgeColor(2,:) );%model
               hold off
               legend( [strcat("model: no treatment", strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup) )))  , ...
                strcat("data 2", " 5FU(i.p. 20mg/kg) are given twice a week for 3 weeks", ",Cell line ",CellLine(2) ),...
                strcat("data 2", "  control group") ,...
                strcat( "Model simulation of data 2 ", strjoin( append ( parameter_label(end-2:end),"=",  string(theta_ControlGroup ), parameter_unit_label(end-2:end), " " )) ) ,...
                ], ...
                'Location', 'best', 'FontSize',12) 

            else %normal process %plot the two cases together 
                for i = 1: Sheets_Dataset_Num %plot treated groups data; One sheet, one dataset
                    Data_20mgkg_TwiceAWeek = cell2mat( Data_20mgkg_TwiceAWeek_cell(i));
                    % plot( Data_20mgkg_TwiceAWeek(:,1), Data_20mgkg_TwiceAWeek(:,2) ,'^','MarkerSize',8  ,'MarkerEdgeColor',MarkerEdgeColor(i,:) );
                    errorbar( Data_20mgkg_TwiceAWeek(:,1), Data_20mgkg_TwiceAWeek(:,2) ,Data_20mgkg_TwiceAWeek(:,4),'^','MarkerSize',8  ,'MarkerEdgeColor',MarkerEdgeColor(i,:) );
                end
                hold on
                for i = 2: Sheets_Dataset_Num %plot control groups data
                    Data_20mgkg_TwiceAWeek = cell2mat( Data_20mgkg_TwiceAWeek_cell(i));
                    plot( Data_20mgkg_TwiceAWeek(:,1), Data_20mgkg_TwiceAWeek(:,3) ,'*','MarkerSize',8,'MarkerEdgeColor',MarkerEdgeColor(i,:) );
                end
                hold on
                plot(Tv_Tumor_data1,  Cfit_Tumor_data1(:,1)/TV_initial_col(1) , 'LineWidth',2, 'color',MarkerEdgeColor(1,:) );%model
                %hold on
                plot(Tv_Tumor_data2,  Cfit_Tumor_variation(:,1)/TV_initial_col(2) , 'LineWidth',2,'color',MarkerEdgeColor(2,:) );%model 
                hold on
                plot(Tv_Tumor_data3,  Cfit_Tumor_data3(:,1)/TV_initial_col(3) , 'LineWidth',2 ,'color',MarkerEdgeColor(3,:) );%model
                hold off           
                legend( [strcat("model: no treatment", strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup) )))  , ...
                append("data ", string(2: Sheets_Dataset_Num), "5FU(i.p. 20mg/kg) are given twice a week for 3 weeks", ",Cell line ",CellLine(2:3)),...
                append("data", string(2: Sheets_Dataset_Num) , "  control group") ,...
                strcat(     'Model simulation of data 2', strjoin( append ( parameter_label(end-2:end),"=",  string( theta_Tv_data2), parameter_unit_label(end-2:end), " " ))),...
                strcat( "Model simulation of data 3" ,strjoin( append ( parameter_label(end-2:end),"=",  string( theta_Tv_data3), parameter_unit_label(end-2:end), " " )) ),...
                ], ...
                'Location', 'best', 'FontSize',12) 
            end 
            hold off 
            grid off;box off
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
%             title("Tumor Volume" +newline+...
%                 "Model parameters are generated by fitting to  data:  " + num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
%                 "Parameter for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " )),'FontSize',12 )
            
            
            
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_elongated(4,:)); %20 mg/kg
            Data_20mgkg_TheOtherDay_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'20mgkg_EveryTheOtherDay.xls');
            Sheets_Dataset_Num = size(Data_20mgkg_TheOtherDay_cell,1);
            CellLine =  ["SW620","HT116"];
            figure('Position', [258.6,93,605.6,470.4]);  
            AllModel_info   = CollectDataIntoNestedStructure ;
            if contains( Seed_Model_DisplayLabel, 'SW620') %only plot SW620 %main 5 to 7
               ErrorBar = AllModel_info(7).ErrorBar ; 
               plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial ,'--', 'LineWidth',2 , 'color', MarkerEdgeColor(3,:) );%model in the absence of treatment
               hold on
               Data_20mgkg_TheOtherDay = cell2mat(Data_20mgkg_TheOtherDay_cell(1));
               errorbar( Data_20mgkg_TheOtherDay(:,1),  Data_20mgkg_TheOtherDay(:,2),ErrorBar(:,1),'^','MarkerEdgeColor',MarkerEdgeColor(1,:) );
               hold on
               errorbar( Data_20mgkg_TheOtherDay(:,1), Data_20mgkg_TheOtherDay(:,3) ,ErrorBar(:,1), '*','MarkerEdgeColor',MarkerEdgeColor(3,:) );
               hold on 
               plot(Tv_Tumor_variation, Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2, 'color',MarkerEdgeColor(1,:) );%model
               hold off
               legend( ["model: no treatment" , ...
                strcat("data- ", "Protocol 7")  %",Cell line ",CellLine(1) ),...
                strcat("data -", "control group") ,...
                  "Model simulation "  ], ...
                'Location', 'best', 'FontSize',12)  
            else
                plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2 );%model
                hold on
                for i = 1: Sheets_Dataset_Num  %plot treated groups
                    Data_20mgkg_TheOtherDay = cell2mat(Data_20mgkg_TheOtherDay_cell(i));
                    plot( Data_20mgkg_TheOtherDay(:,1) , Data_20mgkg_TheOtherDay(:,2) ,'^','MarkerSize',8);
                end
                hold on
                plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial , 'LineWidth',2 );%model in the absence of treatment
                for i = 1: Sheets_Dataset_Num %plot control groups
                    Data_20mgkg_TheOtherDay = cell2mat( Data_20mgkg_TheOtherDay_cell(i));
                    plot( Data_20mgkg_TheOtherDay(:,1) , Data_20mgkg_TheOtherDay(:,3) ,'*','MarkerSize',8);
                end
                hold off
                legend( [  Dose_info_display(4) ,  append("data", string(1: Sheets_Dataset_Num) ,...
                    "5FU(i.p. 20mg/kg) are administered every other day for 2 weeks" ,",Cell line: ", CellLine ),...
                    "model no treatment" ,  append("data", string(1: Sheets_Dataset_Num) , "  control group" ) ], ...
                    'Location', 'best', 'FontSize',12)
            end
            grid off;box off
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
%             title("Tumor Volume" +newline+...
%                 "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
%                 "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
%                 "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup)) ),'FontSize',12 )
%             
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_elongated(7,:)); %20 mg/kg
            Data_20mgkg_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'20mgkg_3timesAweek_3weeks.xls');
            Sheets_Dataset_Num = size( Data_20mgkg_cell,1);
            legend_label = [  Dose_info_display(7) ,...
                append("data", string(1: Sheets_Dataset_Num) ," (i.p. 20mg/kg) are given three times a week for 3 weeks,Cell line HT29") ,...
                "model no treatment" ,  append("data", string(1: Sheets_Dataset_Num) , "  control group" ) ] ;
            figure %2  every other day
            plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2 );%model
            hold on
            for i = 1: Sheets_Dataset_Num  %plot treated groups
                Data_20mgkg = cell2mat(Data_20mgkg_cell(i));
                %plot( Data_20mgkg(:,1) , Data_20mgkg(:,2) ,'^','MarkerSize',8);
                errorbar(Data_20mgkg(:,1) , Data_20mgkg(:,2), Data_20mgkg(:,4) ,'Marker', '^','MarkerSize',8                )
            end
            hold on
            plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial , 'LineWidth',2 );%model in the absence of treatment
            for i = 1: Sheets_Dataset_Num %plot control groups
                Data_20mgkg = cell2mat( Data_20mgkg_cell(i));
                %plot( Data_20mgkg(:,1) , Data_20mgkg(:,3) ,'*','MarkerSize',8);
                errorbar(Data_20mgkg(:,1) , Data_20mgkg(:,3),  Data_20mgkg(:,5) ,'Marker', '^','MarkerSize',8                )
                hold on
            end
%             if FeaturedDose == 20
%                 origin_para = [   580.589661618  0.00445 0.24860 0.26015  0.000282501  1.883650695 0.00006 ];
%                 [Tv_Tumor_OrigPmax,Cfit_Tumor_OrigPmax]  = myfun_kinetics( origin_para ,TV_initial, Dose_info_elongated(7,:)); %20 mg/kg
%                 plot(Tv_Tumor_OrigPmax, Cfit_Tumor_OrigPmax(:,1)/TV_initial  ,'-.' ,'LineWidth',2 )
%                 legend_label(end + 1) = "simulaion using model 4";
%             end
            hold off
            grid
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
            
            legend(  legend_label,...
                'Location', 'best', 'FontSize',12)
            title("Tumor Volume" +newline+...
                "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
                "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
                "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup)) ),'FontSize',12 )
            
            
            
            
            
            
            
            
        elseif FocusOnTumorVolume_SeparatePlots_CurrentDose == 100
            Data_100mgkg_Weekly_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'TV_100mgkgweekly.xls');
            Data_100mgkg_Weekly = cell2mat( Data_100mgkg_Weekly_cell);
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_elongated(end-1,:));
            figure
            plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2 );%model
            hold on
            plot( Data_100mgkg_Weekly (:,1) , Data_100mgkg_Weekly (:,2) ,'^','MarkerSize',8);
            
            hold on
            plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial , 'LineWidth',2 );%model in the absence of treatment
            hold on
            plot( Data_100mgkg_Weekly (:,1) ,  Data_100mgkg_Weekly (:,3) ,'*','MarkerSize',8);
            
            hold off
            grid
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
            legend( [  Dose_info_display(end-1) ,   "data 1"+"  5FU(i.p. 100 mg/kg) are administered weekly for four times, Cell line: Colon 38",...
                "model: no treatment" ,  "data 1: control group"    ], ...
                'Location', 'best', 'FontSize',12)
            title("Tumor Volume" +newline+...
                "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
                "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
                "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup)) ),'FontSize',12 )
            
        elseif FocusOnTumorVolume_SeparatePlots_CurrentDose == 50
            Data_50mgkg_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'50mgkg_daily_percentage.xls');
            Data_50mgkg_array = cell2mat( Data_50mgkg_cell);
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_elongated(6,:));
            figure
            plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2 );%model
            hold on
            plot( Data_50mgkg_array (:,1) , Data_50mgkg_array(:,2) ,'^','MarkerSize',8);
            hold on
            plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial , 'LineWidth',2 );%model in the absence of treatment
            hold on
            plot( Data_50mgkg_array(:,1) ,  Data_50mgkg_array(:,3) ,'*','MarkerSize',8);
            hold off
            grid
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
            legend( [  Dose_info_display(6) ,   "data 1"+"  5FU(i.p. 50 mg/kg) are administered daily for 4 times, Cell line: Colon 26",...
                "model: no treatment" ,  "data 1: control group"    ], ...
                'Location', 'best', 'FontSize',12)
            title("Tumor Volume" +newline+...
                "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
                "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
                "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup)) ),'FontSize',12 )
            
            
            
            
            
        elseif FocusOnTumorVolume_SeparatePlots_CurrentDose == 10
            Data_10mgkg_3TimesAweek_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'10mgkg_3timesAweek_3weeks_LabelWeek.xls');
            Data_10mgkg_2TimesAweek_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'10mgkg_TwiceAweek_5week_LabelWeek.xls');
            Sheets_Dataset_Num_1 = size(Data_10mgkg_3TimesAweek_cell,1);
            Sheets_Dataset_Num_2 = size(Data_10mgkg_2TimesAweek_cell,1);
            
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_elongated(2,:));
            figure
            plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2 );%model
            hold on
            for i = 1: Sheets_Dataset_Num_1  %plot treated groups
                Data_10mgkg_3TimesAweek = cell2mat(Data_10mgkg_3TimesAweek_cell(i));
                plot( Data_10mgkg_3TimesAweek(:,1) , Data_10mgkg_3TimesAweek(:,2) ,'^','MarkerSize',8);
            end
            hold on
            plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial , 'LineWidth',2 );%model in the absence of treatment
            hold on
            for i = 1: Sheets_Dataset_Num_1 %plot control groups
                Data_10mgkg_3TimesAweek = cell2mat( Data_10mgkg_3TimesAweek_cell(i));
                plot(  Data_10mgkg_3TimesAweek(:,1) ,  Data_10mgkg_3TimesAweek(:,3) ,'*','MarkerSize',8);
            end
            hold off
            grid
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
            legend( [  Dose_info_display(2) ,   "data 1"+"  5FU(i.p. 10 mg/kg) are administered three times a week for 3 weeks, Cell line:HT29",...
                "model: no treatment" ,  "data 1: control group"    ], ...
                'Location', 'best', 'FontSize',12)
            title("Tumor Volume" +newline+...
                "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
                "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
                "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup)) ),'FontSize',12 )
            
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_elongated(3,:));
            figure
            plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2 );%model
            hold on
            for i = 1: Sheets_Dataset_Num_2  %plot treated groups
                Data_10mgkg_weekly = cell2mat(Data_10mgkg_2TimesAweek_cell (i));
                plot( Data_10mgkg_weekly(:,1) , Data_10mgkg_weekly(:,2) ,'^','MarkerSize',8);
            end
            hold on
            plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial , 'LineWidth',2 );%model in the absence of treatment
            for i = 1: Sheets_Dataset_Num_2 %plot control groups
                Data_10mgkg_weekly = cell2mat(Data_10mgkg_2TimesAweek_cell (i));
                plot( Data_10mgkg_weekly(:,1) , Data_10mgkg_weekly(:,3) ,'*','MarkerSize',8);
            end
            hold off
            grid
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
            legend( [ Dose_info_display(3), "data 1: 5FU(i.p. 10 mg/kg) are administered twice a week for 5 weeks (nonspecific cell line)",...
                "model no treatment" ,  "data 1: control group"], ...
                'Location', 'best', 'FontSize',12)
            title("Tumor Volume" +newline+...
                "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
                "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
                "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup) )) , 'FontSize',12 )
            
            
        elseif FocusOnTumorVolume_SeparatePlots_CurrentDose == 5
            [Tv_Tumor_variation,Cfit_Tumor_variation]  = myfun_kinetics(theta_Tv,TV_initial, Dose_info_elongated(1,:));
            Data_5mgkg_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'5mgkg_daily.xls');
            Treatement_error = [0.2150    0.2087    0.3665    0.4910    0.3751    0.3313    0.7134    1.0437    0.62];
            control_error = [0.3896    0.2980    0.3335    0.4938    0.4131    0.8720    1.2900    2.5168    2.7934];
            
            figure('Position', [258.6,93,605.6,470.4]);  %1  every the other week
            plot(Tv_Tumor_variation,  Cfit_Tumor_variation(:,1)/TV_initial , 'LineWidth',2 );%model
            hold on
            Data_5mgkg_array = cell2mat(Data_5mgkg_cell(1));
            errorbar( Data_5mgkg_array(:,1) , Data_5mgkg_array(:,2) , Treatement_error ,'o', 'MarkerSize',8);
            hold on
            plot(Tfit_Absolute_ControlGroup, Cfit_Absolute_ControlGroup(:,1)/TV_initial , 'LineWidth',2 );%model in the absence of treatment
            hold on
            errorbar( Data_5mgkg_array(:,1) , Data_5mgkg_array(:,3) , control_error ,'^', 'MarkerSize',8); 
            hold off
            grid off;box off
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume(fold change)')
            xticks( 0: Dose_info_elongated(end,3): fix(   Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )*Dose_info_elongated(end,3) )
            xticklabels(  num2cell( 0 :fix(  Dose_info_elongated(end,1)/Dose_info_elongated(end,3) )  ) )
            legend( [  "Model simulation",  "data 1" + "Protocol 5",...
                "model: no treatment" ,  "data 1: control group" ], ...
                'Location', 'best', 'FontSize',12)
%             title("Tumor Volume" +newline+...
%                 "Model parameters are generated by fitting to " +" data: "+ num2str(FeaturedDose)+" mg/kg "+ Seed_Model_DisplayLabel+newline+...
%                 "Parameter estimates for treated group: " +strjoin( append ( parameter_label,"=",  string( theta_Tv), parameter_unit_label, " " ))+newline+...
%                 "Parameters estimates for control group: " +strjoin( append ( parameter_label_Control ,"=",  string( theta_ControlGroup)) ),'FontSize',12 )
        end
    end
elseif FocusOnTumorVolume == 0 %plot intermediates
    [Dose_info_uniq,~,~] = unique(Dose_info,'rows','stable');
    if   IntermediateSubplot_InOneFig == 1
        Elongation = 5*day2min;
        Dose_info_uniq(:,1) = Dose_info_uniq(:,1)+ Elongation ;
        for i= 1: size(  Dose_info_uniq,1  )
            Plot_IntermediateSubplot(Dose_info_uniq(i,:) , Dose_info_display(i) );
        end
    elseif   IntermediateSubplot_InOneFig == 0
        Dose_info_uniq(:,1) = (Dose_info_uniq(end,1)+ Elongation )*ones(size(  Dose_info_uniq(:,1)  ) );
        NumCases = size(Dose_info_uniq,1);
        %One cell for the realization of one row of Dose_info
        %Set apart time points from the vlaues of variables
        T_up2DSB_col = cell( NumCases,1   );
        C_up2DSB_col = cell( NumCases,1   );
        for  i = 1: NumCases-1%5, 10, 20,50, 100
            [T_up2DSB, C_up2DSB] = kinetics_uptoDSB( Dose_info_uniq(i,:) );
            T_up2DSB_col{i} = T_up2DSB;
            C_up2DSB_col{i} = C_up2DSB;
        end
        shape = ['o' ,'+', '*'];
        t_Interstitial = cell2mat(T_allsheets_output_cell(1) ) ;
        c_Interstitial = cell2mat(C_allsheets_output_cell(1) ) ;
        figure
        for i = 1: NumCases-1   %model , except for control
            T_up2DSB = cell2mat( T_up2DSB_col(i) );
            C_up2DSB = cell2mat( C_up2DSB_col(i) );
            if TruncateDuration_ToDataTimeFrame ==1
                [~,idx] = min(abs( T_up2DSB- max(t_Interstitial)-100 ));%only one
            else
                [~,idx] = min(abs( T_up2DSB- 2*week2min ));%only one
            end
            plot( T_up2DSB(1:idx), C_up2DSB(1:idx,1)/UnitConversion , 'LineWidth',2,'DisplayName',Dose_info_display(i) ); %model
            hold on
        end
        for i = 1: size( dose_c0,1 ) -1   %plot central Interstitial data of different doses(5, 10, 20),except for 50
            data_plt =  plot(t_Interstitial_5and10and20(:,i), c_Interstitial_5and10and20(:,i)  ,shape(i) , 'MarkerSize',8, 'DisplayName', strcat( 'data' , DoseLabel(i),' one dose' ));%data
            hold on
        end
        data_plt =  plot(t_Interstitial, c_Interstitial  , '^', 'MarkerSize',8, 'DisplayName', strcat( 'data' , '100mg/kg one dose' ));%data for 100
        hold off
        grid
        xlabel('Time(min)')
        ylabel('concentration(μg/mL)')
        title('Central Interstitial')
        legend('Location', 'best','FontSize',11)
        
        species = 5; %5FU anabolites
        figure
        for i = 1: NumCases-1   %model
            T_up2DSB = cell2mat( T_up2DSB_col(i) );
            C_up2DSB = cell2mat( C_up2DSB_col(i) );
            if TruncateDuration_ToDataTimeFrame ==1
                idx_col= find( C_up2DSB(:,species) == min(  C_up2DSB(:,species)   )  ); %min should be 0; more than one
                idx =  idx_col(2);
            else
                [~,idx] = min(abs( T_up2DSB- 2*week2min ));
            end
            plot( T_up2DSB(1: idx), C_up2DSB(1: idx,species) , 'LineWidth',2,'DisplayName', Dose_info_display(i)); %model
            hold on
        end
        hold off
        grid
        xlabel('Time(h)')
        xticks( 0: 10*60: T_up2DSB( idx) )
        xticklabels(  num2cell(0: 10: fix(T_up2DSB( idx)/60  ) ) )
        ylabel('amount(pmol/mg tissue)')
        title('5FU anabolites')
        legend('Location', 'best','FontSize',11)
        
        t_genome = cell2mat(T_allsheets_output_cell(4) ) ;
        species = 6; %F-RNA
        c_RNA  = cell2mat(C_allsheets_output_cell(4));
        figure
        for i = 1: NumCases-1 %model
            T_up2DSB = cell2mat( T_up2DSB_col(i) );
            C_up2DSB = cell2mat( C_up2DSB_col(i) );
            if TruncateDuration_ToDataTimeFrame ==1
                [~,idx] = min(abs( T_up2DSB- 14*day2min)); %only one
            else
                [~,idx] = min(abs( T_up2DSB- 2*week2min ));
            end
            plot( T_up2DSB(1: idx), C_up2DSB(1:idx,species ) , 'LineWidth',2,'DisplayName', Dose_info_display(i)); %model
            hold on
        end
        plot( t_genome, c_RNA,'^', 'MarkerSize',7,'MarkerFaceColor', [0.8500,0.3250,0.0980],'DisplayName','data 100 mg/kg');%data
        hold off
        grid
        xlabel('Time(day)')
        xticks([1;2;3;7;10; 14 ]*day2min)
        xticklabels({'1','2','3','7','10','14'} )
        ylabel('amount(pmol/mg tissue)')
        title('F-RNA')
        legend('Location', 'best','FontSize',11)
        
        species = 7; %F-DNA
        c_DNA  = cell2mat(C_allsheets_output_cell(5));
        figure
        for i = 1:NumCases-1 %model
            T_up2DSB = cell2mat( T_up2DSB_col(i) );
            C_up2DSB = cell2mat( C_up2DSB_col(i) );
            if TruncateDuration_ToDataTimeFrame ==1
                [~,idx] = min(abs( T_up2DSB- 14*day2min));
            else
                [~,idx] = min(abs( T_up2DSB- 2*week2min ));
            end
            plot( T_up2DSB(1: idx), C_up2DSB(1:idx,species ) , 'LineWidth',2,'DisplayName', Dose_info_display(i)); %model
            hold on
        end
        plot( t_genome,  c_DNA,'^', 'MarkerSize',7,'MarkerFaceColor', [0.8500,0.3250,0.0980],'DisplayName','data 100 mg/kg');%data
        hold off
        grid
        xlabel('Time(day)')
        xticks([1;2;3;7;10; 14 ]*day2min)
        xticklabels({'1','2','3','7','10','14'} )
        ylabel('amount(pmol/mg tissue)')
        title('F-DNA')
        legend('Location', 'best','FontSize',11)
        
        t_TS = cell2mat(T_allsheets_output_cell(6) ) ;
        [T_data_xtick, T_data_xlabel] = XtickLabelCreation(t_TS  ,t_TS /60 ); %T_data , T_data/24/60
        figure
        for i = 1:NumCases-1 %model
            T_up2DSB = cell2mat( T_up2DSB_col(i) );
            C_up2DSB = cell2mat( C_up2DSB_col(i) );
            if TruncateDuration_ToDataTimeFrame ==1
                [~,idx] = min(abs( T_up2DSB- t_TS(end) -10*60  ));
            else
                [~,idx] = min(abs( T_up2DSB- 2*week2min ));
            end
            plot( T_up2DSB(1: idx), C_up2DSB(1:idx,end )*100 , 'LineWidth',2,'DisplayName', Dose_info_display(i)); %model
            hold on
        end
        hold off
        grid
        xlabel('Time(h)')
        xticks(T_data_xtick)
        xticklabels(T_data_xlabel)
        ylabel('Percentage')
        title(' %free TS')
        legend('Location', 'best','FontSize',11)
        
        t_dNTP = cell2mat(T_allsheets_output_cell(9) ) ;
        [T_data_xtick, T_data_xlabel] = XtickLabelCreation( t_dNTP , t_dNTP/60 ); %T_data , T_data/24/60
        species = 10;
        figure
        for i = 1:NumCases-1 %model
            T_up2DSB = cell2mat( T_up2DSB_col(i) );
            C_up2DSB = cell2mat( C_up2DSB_col(i) );
            if TruncateDuration_ToDataTimeFrame ==1
                [~,idx] = min(abs( T_up2DSB- t_dNTP(end)-10*60  ));
            else
                [~,idx] = min(abs( T_up2DSB- 2*week2min ));
            end
            plot( T_up2DSB(1: idx), C_up2DSB(1:idx,species  ) , 'LineWidth',2,'DisplayName', Dose_info_display(i)); %model
            hold on
        end
        hold off
        grid
        xlabel('Time(h)')
        xticks(T_data_xtick)
        xticklabels(T_data_xlabel)
        ylabel('Perturbation')
        title(' dNTP pool imabalance')
        legend('Location', 'best','FontSize',11)
        
        t_DSB = cell2mat(T_allsheets_output_cell(10) ) ;
        [T_data_xtick, T_data_xlabel] = XtickLabelCreation( t_DSB , t_DSB/60 ); %T_data , T_data/24/60
        species = 11;
        figure
        for i = 1:NumCases-1 %model
            T_up2DSB = cell2mat( T_up2DSB_col(i) );
            C_up2DSB = cell2mat( C_up2DSB_col(i) );
            if TruncateDuration_ToDataTimeFrame ==1
                [~,idx] = min(abs( T_up2DSB- t_DSB(end)-10*60  ));
            else
                [~,idx] = min(abs( T_up2DSB- 2*week2min ));
            end
            plot( T_up2DSB(1: idx), C_up2DSB(1:idx,species  ) , 'LineWidth',2,'DisplayName', Dose_info_display(i)); %model
            hold on
        end
        hold off
        grid
        xlabel('Time(h)')
        xticks(T_data_xtick)
        xticklabels(T_data_xlabel)
        ylabel('γH2AX foci count(thousands)')
        title('DSB')
        legend('Location', 'best','FontSize',11)
    end
end
end

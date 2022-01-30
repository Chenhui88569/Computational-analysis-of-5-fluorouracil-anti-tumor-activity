%100mg/kg VS 20mg/kg
%100 weekly for 2 weeks VS 20mg/kg Twice a week for 5 weeks
%Function responsible for plotting @ Tv_plot
%(Total dose should be the same )
%the length of therapy should be the same.

% The purpose of comparison or the useful information that can give us
% from this study is that multiple doses of low amount of drug can exhibit
% more cytotoxic effects or generate more surpressive effects on tumor
% growth as compared to the single dose of higher amount of drug. The
% precondition is that the length of treatment for the two cases is the
% same. 
clear all
close all
UnitConversion = [  10^6/130.077, 24*60,10^-3 ,7*24*60   ];
UnitConversion_cell = num2cell(UnitConversion  );
[Plasma_UnitConversion, day2min  ,mm3tocm3  ,week2min  ] = UnitConversion_cell{:};
Elongation  =  20*day2min;
SameSchedule = 0;
    if  SameSchedule == 0
        DurInweek_100 = 1;
        DurInweek_20 = 1;
        DurInweek_50 = 1;
        FreqEveryWeek_100= 1;
        FreqEveryWeek_20 = 5;
        FreqEveryWeek_50 = 2;
    elseif  SameSchedule == 1 %haven'e been used so far
        DurInweek_100 = 2;
        DurInweek_20 = 2;
        FreqEveryWeek_100= 2;
        FreqEveryWeek_20 = 2;
    end
    Seed_Schedule_DisplayLabel = [ strcat("100 mg/kg " , num2str( FreqEveryWeek_100), " time(s)",  "a week for ",  num2str(DurInweek_100), " week(s)" ),...
        strcat( "20mg/kg ", num2str( FreqEveryWeek_20), " time(s)",  "a week for ", num2str(DurInweek_20) ," week(s)"),...
         strcat( "50mg/kg ", num2str( FreqEveryWeek_50), " time(s)",  "a week for ", num2str(DurInweek_50) ," week(s)") ]; %string array
    Seed_Schedule_pool = [DurInweek_100*week2min   DurInweek_100*FreqEveryWeek_100-1	 7*day2min/FreqEveryWeek_100   92.9936306*Plasma_UnitConversion
        DurInweek_20*week2min	 DurInweek_20*FreqEveryWeek_20-1	7*day2min/FreqEveryWeek_20   93/5*Plasma_UnitConversion
                DurInweek_50*week2min	 DurInweek_50*FreqEveryWeek_50-1	7*day2min/FreqEveryWeek_50   93/2*Plasma_UnitConversion ];
 
    Total_dose = [100;20;50].* (Seed_Schedule_pool(:,2)+1);
    Seed_Schedule_pool_elongation = Seed_Schedule_pool;
    Seed_Schedule_pool_elongation(:,1) = (Seed_Schedule_pool_elongation(end,1)+ Elongation )*ones(size( Seed_Schedule_pool_elongation(:,1)  ) ); %Duration are the same
    
    Data_cell = readtable('.\TV_100mgkgweekly.xls');
    Data_array = table2array(Data_cell);
    Data_array(:,1)  = Data_array(:,1)*day2min;
    Control  = 125*mm3tocm3; % mm^3
    TV_initial  = Data_array(1,2)*Control ;
    %ParaFromTreatedGroupData_Dose : 100
    theta_TV_ControlGroup_100 = [3.3616  9e-05 0.35014];
    %["T_{tumor}","E_{max}","EC_{50}", "\lambda_{0}","\gamma","k_{out}","\gamma_{hill}"];
    %theta_Tv_treated_100 = [ 1124.73676;36.55462;1515.20249 ;0.00030;1.56684;0.00002 ];%without sigmoidity parameter
    theta_Tv_treated_100 =  [1407.45420;136.49992;12.20191;0.00042;1.77324;0.00005];
    theta_Tv_treated_100_cell  = num2cell(theta_Tv_treated_100 );
    [lag_tumor_100, E_max_100, EC50_100, lambda_0_100 , gamma_100,k_out_100   ] = theta_Tv_treated_100_cell{:};
    theta_Tv_treated_100(end+1) = theta_TV_ControlGroup_100(1);
    theta_Tv_treated_100(end+1)  = 3;
    
    theta_Tv_treated_20 =   [1696.12850; 22.23930;8.26760;0.00020;0.02914;0.00010 ];  
      
    theta_Tv_treated_20_cell  = num2cell(theta_Tv_treated_20 );
    [lag_tumor_20, E_max_20, EC50_20, lambda_0_20 , gamma_20,k_out_20   ] = theta_Tv_treated_20_cell{:};
    theta_TV_ControlGroup_20 = [1.4466  0.01848  9.1215]; % ["\chi_{max}", "\lambda_{0}","\gamma"];\
    theta_Tv_treated_20 (end+1) = theta_TV_ControlGroup_20(1);
    theta_Tv_treated_20 (end+1) =3; %gamma hill
    
    theta_Tv_treated = [1490.36176;60.58025;48.16340;0.00009;0.22718;0.00003];
    theta_Tv_treated_cell  = num2cell(theta_Tv_treated );
    [lag_tumor, E_max, EC50, lambda_0 , gamma,k_out  ] = theta_Tv_treated_cell{:}; 
    theta_TV_ControlGroup = [3.55035  0.00005 0.26759]; % ["\chi_{max}", "\lambda_{0}","\gamma"];\
    theta_Tv_treated(end+1) = theta_TV_ControlGroup(1);
    theta_Tv_treated(end+1) =3; %gamma hill

    %Seed_schedule used in comparision, plot tumor volume according to the
    %schedule first
    myfun_kinetics  = @kinetics_TvPlot_DoseDifference; %kinetics_TvPlot_DoseDifference(theta,TV_initial, Dose_info)
    [Tfit_100,Cfit_100]  = myfun_kinetics(theta_Tv_treated_100,TV_initial, Seed_Schedule_pool_elongation(1,:));
    [Tfit_20,Cfit_20]  = myfun_kinetics(theta_Tv_treated_20,TV_initial,Seed_Schedule_pool_elongation(2,:)); %20 mg/kg
    [Tfit_50,Cfit_50]  = myfun_kinetics(theta_Tv_treated,TV_initial, Seed_Schedule_pool_elongation(3,:));
 
    figure
    plot( Tfit_100,  Cfit_100/TV_initial ,'LineWidth',2)
    hold on
    plot( Tfit_20,  Cfit_20/TV_initial ,'LineWidth',2    )
    hold off
    legend(Seed_Schedule_DisplayLabel ,'FontSize',12)
    grid on
    
    [T_plasma_100, C_uptoDSB_100] = kinetics_uptoDSB( Seed_Schedule_pool(1,:) );
    T_plasma_100_day =T_plasma_100/day2min; 
    C_plasma_100 = C_uptoDSB_100(:,1)/Plasma_UnitConversion;
    [T_plasma_20, C_uptoDSB_20] = kinetics_uptoDSB( Seed_Schedule_pool(2,:) );
    T_plasma_20_day = T_plasma_20/day2min;
    C_plasma_20 = C_uptoDSB_20(:,1)/Plasma_UnitConversion;
    [T_plasma_50, C_uptoDSB_50] = kinetics_uptoDSB( Seed_Schedule_pool(3,:) );
    T_plasma_50_day = T_plasma_50/day2min;
    C_plasma_50 = C_uptoDSB_50(:,1)/Plasma_UnitConversion;

    
    MIC = 5;
    AUC_plasma_100 = trapz( T_plasma_100_day   , C_plasma_100 ); %Trapezoidal numerical integration
    AUC_plasma_100_MIC = AUC_plasma_100/MIC  ;
    AUC_plasma_20 = trapz(T_plasma_20_day    ,C_plasma_20  );
    AUC_plasma_20_MIC = AUC_plasma_20/MIC  ;
    AUC_plasma_50 = trapz(T_plasma_50_day    ,C_plasma_50  );
    AUC_plasma_50_MIC = AUC_plasma_50/MIC  ;

    
    Find_MIC_idx_100  =  find( C_plasma_100>MIC );
    MIC_idx_100_last = Find_MIC_idx_100(end);
    Find_MIC_idx_20  =  find( C_plasma_20>MIC );
    MIC_idx_20_last = Find_MIC_idx_20(end); 
        Find_MIC_idx_50  =  find( C_plasma_50>MIC );
    MIC_idx_50_last = Find_MIC_idx_50(end); 
    TimeAboveMIC_fun = @(T_plasma ,C_plasma,idx)  interp1(C_plasma(idx:idx+1),T_plasma(idx:idx+1) , MIC   ); 
    TimeAboveMIC_100 = TimeAboveMIC_fun(T_plasma_100_day ,C_plasma_100, MIC_idx_100_last  );
    TimeAboveMIC_20 = TimeAboveMIC_fun(T_plasma_20_day ,C_plasma_20, MIC_idx_20_last  );
    TimeAboveMIC_50 = TimeAboveMIC_fun(T_plasma_50_day ,C_plasma_50, MIC_idx_50_last  );

    E_drug_100 = @(c_plasma) E_max_100.*c_plasma.^3./( EC50_100.^3 +c_plasma.^3   ); %obj has imaginary part
    E_drug_20 = @(c_plasma) E_max_20*c_plasma.^3 ./( EC50_20.^3+c_plasma.^3  ); %obj has imaginary part
    E_drug_50 = @(c_plasma) E_max*c_plasma.^3 ./( EC50.^3+c_plasma.^3  ); %obj has imaginary part

    E_drug = @(c_plasma) E_max.*c_plasma.^3./( EC50.^3 +c_plasma.^3   ); %obj has imaginary part
    disp( strcat( "Suppose dose regimen are ",   Seed_Schedule_DisplayLabel(1) ," and " , Seed_Schedule_DisplayLabel(2) ) );
    fprintf('\t Total dose for 100 mg/kg  = %8.5f\n',Total_dose(1));
    fprintf('\t Total dose for 20 mg/kg  = %8.5f\n',Total_dose(2));    
    fprintf('\t Total dose for 50 mg/kg  = %8.5f\n',Total_dose(3));

    fprintf('\t Time>MIC for 100 mg/kg  = %8.5f\n',TimeAboveMIC_100);
    fprintf('\t Time>MIC for 20 mg/kg  = %8.5f\n',TimeAboveMIC_20);
    fprintf('\t Time>MIC for 50 mg/kg  = %8.5f\n',TimeAboveMIC_50);
   
    fprintf("\t Area under the blood (or plasma) concentration–time curve for "+ Seed_Schedule_DisplayLabel(1) +"  = %8.5f\n",AUC_plasma_100 );
    fprintf(" \t Area under the blood (or plasma) concentration–time curve for " + Seed_Schedule_DisplayLabel(2) +"  = %8.5f\n",AUC_plasma_20  );
    fprintf(" \t Area under the blood (or plasma) concentration–time curve for " + Seed_Schedule_DisplayLabel(3) +"  = %8.5f\n",AUC_plasma_50  );
    disp(" 1. If model parameters that are based on 100 mg/kg dosage are used")
   
    AUC_TV_100 = trapz(T_plasma_100_day   , E_drug_100( C_plasma_100   )); %Trapezoidal numerical integration
    AUC_TV_20 = trapz(T_plasma_20_day   ,E_drug_100( C_plasma_20   ));
    AUC_TV_50 = trapz(T_plasma_50_day   ,E_drug_100( C_plasma_50   ));
    fprintf("\t\t  Area under the E_max curve for "+ Seed_Schedule_DisplayLabel(1) +"  = %8.5f\n",AUC_TV_100 );
    fprintf(" \t\t Area under the E_max curve for " + Seed_Schedule_DisplayLabel(2) +"  = %8.5f\n",AUC_TV_20  );
    fprintf(" \t\t Area under the E_max curve for " + Seed_Schedule_DisplayLabel(3) +"  = %8.5f\n",AUC_TV_50  );
  
    disp(" 2. If model parameters that are based on 20 mg/kg dosage are used")
    AUC_TV_100 = trapz(T_plasma_100_day   ,E_drug_20( C_plasma_100   )); %Trapezoidal numerical integration
    AUC_TV_20 = trapz(T_plasma_20_day   ,E_drug_20( C_plasma_20   ));
    AUC_TV_50 = trapz(T_plasma_50_day   ,E_drug_20( C_plasma_50   ));

    fprintf("\t\t Area under the E_max curve for "+ Seed_Schedule_DisplayLabel(1) +"  = %8.5f\n",AUC_TV_100 );
    fprintf(" \t\t Area under the E_max curve for " + Seed_Schedule_DisplayLabel(2) +"  = %8.5f\n",AUC_TV_20  );
    fprintf(" \t\t Area under the E_max curve for " + Seed_Schedule_DisplayLabel(3) +"  = %8.5f\n",AUC_TV_50  );
    disp(" 3. If model parameters that are based on 100mg/kg and 20 mg/kg dosage respectively are used")
    AUC_TV_100 = trapz(T_plasma_100_day    ,E_drug_100( C_plasma_100   )); %Trapezoidal numerical integration
    AUC_TV_20 = trapz(T_plasma_20_day ,E_drug_20( C_plasma_20   ));
    AUC_TV_50 = trapz(T_plasma_50_day ,E_drug_50( C_plasma_50   ));

    fprintf("\t\t Area under the E_max curve for "+ Seed_Schedule_DisplayLabel(1) +"  = %8.5f\n",AUC_TV_100 );
    fprintf(" \t\tArea under the E_max curve for " + Seed_Schedule_DisplayLabel(2) +"  = %8.5f\n",AUC_TV_20  );
    fprintf(" \t\tArea under the E_max curve for " + Seed_Schedule_DisplayLabel(3) +"  = %8.5f\n",AUC_TV_50  );
    %concentration-effect curve (100,20)
    %NarrowConcentrationTo =60;
    NarrowConcentrationTo = 40;
    figure
    [~,idx_0] = min(abs(C_plasma_100- 0  ));
    [~,idx_max] = min(abs(C_plasma_100(1:idx_0,1)-  NarrowConcentrationTo ));
    concentration_sample = 1:1:120;
    plot(concentration_sample  ,  E_drug_100( concentration_sample   ) ,'LineWidth',2)
    hold on
    [~,idx_0] = min(abs(C_plasma_20- 0  ));
    [~,idx_max] = min(abs(C_plasma_20(1:idx_0,1) - NarrowConcentrationTo   ));
    plot(  concentration_sample, E_drug_20( concentration_sample) , 'LineWidth',2    )
    hold on 
    [~,idx_0] = min(abs(C_plasma_50- 0  ));
    [~,idx_max] = min(abs(C_plasma_50(1:idx_0,1) - NarrowConcentrationTo   ));
    plot(concentration_sample , E_drug_50( concentration_sample ), 'LineWidth',2    )
    hold off
    ylabel(" effect(units) ")
    xlabel("\mu g/mL  ")
    legend(["one dose 100 mg/kg" , "one dose 20 mg/kg", "one dose 50 mg/kg"] ,'FontSize',10,'Location','best')
    title("5FU concentration-effect curve (100 mg/kg, 20mg/kg,50mg/kg)")
    grid on
    
     %time-concentration curve (100, 20) over 1 weeks
    NarrowTimeFrameTo = 2;
    if NarrowTimeFrameTo> DurInweek_100
        NarrowTimeFrameTo  = DurInweek_100;
    end
    [~,idx_100] = min(abs( T_plasma_100 - NarrowTimeFrameTo*week2min ));
    [~,idx_20] = min(abs( T_plasma_20 - NarrowTimeFrameTo*week2min ));
    [~,idx_50] = min(abs( T_plasma_50 - NarrowTimeFrameTo*week2min ));

    figure
    plot(  T_plasma_100(1:idx_100),   C_plasma_100(1:idx_100) ,'LineWidth',2)
    hold on
    plot( T_plasma_20(1:idx_20) ,  C_plasma_20(1:idx_20), 'LineWidth',2    ) 
    hold on 
    plot( T_plasma_50(1:idx_50) ,  C_plasma_50(1:idx_50), 'LineWidth',2    )
    hold off 
    xticks( 0: day2min :  NarrowTimeFrameTo*week2min ) % in day
    xticklabels(  num2cell( 0 :1: NarrowTimeFrameTo*7 )  )
    ylabel(" \mug/mL ")
    xlabel("days")
    legend(Seed_Schedule_DisplayLabel ,'FontSize',10,'Location','best')
    title("5FU time-concentration (20mg/kg, 100mg/kg,50mg/kg)")
    grid on
    
    %time-effect curve (100, 20) over 2 weeks
    figure
    plot(  T_plasma_100 ,  E_drug_100(  C_plasma_100) ,'LineWidth',2)
    hold on
    plot( T_plasma_20 ,  E_drug_20(  C_plasma_20), 'LineWidth',2    )
    hold on 
    plot( T_plasma_50 ,  E_drug_50(  C_plasma_50), 'LineWidth',2    )

    hold off
    xticks( 0: week2min*0.5:  DurInweek_100*week2min )
    xticklabels(  num2cell( 0 :1/2:DurInweek_100 )  )
    ylabel(" effect(units) ")
    xlabel("weeks")
    legend(Seed_Schedule_DisplayLabel ,'FontSize',10,'Location','best')
    title("5FU time-effect curve (20mg/kg, 100mg/kg)")
    grid on
    
    %time-effect curve and time-plasma curve(100)
    NarrowTimeFrameTo = 1;
    if NarrowTimeFrameTo> DurInweek_100
        NarrowTimeFrameTo  =0;
    end
    [~,idx_100] = min(abs( T_plasma_100 - NarrowTimeFrameTo*day2min ));
    figure
    yyaxis left
    plot(  T_plasma_100(1:idx_100) ,  E_drug_100(  C_plasma_100(1:idx_100)) ,'LineWidth',2)
    ylabel("effects(units)")
    hold on
    yyaxis right
    plot( T_plasma_100(1:idx_100) , C_plasma_100(1:idx_100), 'LineWidth',2    )
    ylabel("Plasma concentration(\mug/mL) ")
    hold off
    xticks( 0: day2min*0.5:  NarrowTimeFrameTo*day2min  )
    xticklabels(  num2cell( 0 :1/2:NarrowTimeFrameTo )  )
    xlabel("Time(days)")
    legend(["drug time-effect 100 mg/kg","plasma concentration 100 mg/kg"],'FontSize',10,'Location','best')
    title("5FU time-effect and time-concentration curve(100 mg/kg) ")
    grid on
    
    %time-effect curve and time-plasma curve(20)
    NarrowTimeFrameTo = 1;
    if NarrowTimeFrameTo> DurInweek_100
        NarrowTimeFrameTo  =0;
    end
    [~,idx_100] = min(abs( T_plasma_20 - NarrowTimeFrameTo*day2min ));
    figure
    yyaxis left
    plot(  T_plasma_20(1:idx_100) ,  E_drug_20(  C_plasma_20(1:idx_100)) ,'LineWidth',2)
    ylabel("effects(units)")
    hold on
    yyaxis right
    plot( T_plasma_20(1:idx_100) , C_plasma_20(1:idx_100), 'LineWidth',2    )
    ylabel("Plasma concentration(\mug/mL) ")
    hold off
    xticks( 0: day2min*0.5:  NarrowTimeFrameTo*day2min  )
    xticklabels(  num2cell( 0 :1/2:NarrowTimeFrameTo )  )
    xlabel("Time(days)")
    legend(["drug time-effect 20 mg/kg","plasma concentration 20 mg/kg"],'FontSize',10,'Location','best')
    title("5FU time-effect and time-concentration curve (20 mg/kg)")
    grid on
    
    
    %time-effect curve and time-plasma curve(100,20)
    NarrowTimeFrameTo = 6/24; 
    if NarrowTimeFrameTo> DurInweek_100
        NarrowTimeFrameTo  =0;
    end
    [~,idx_100] = min(abs( T_plasma_100 - NarrowTimeFrameTo*day2min ));
    [~,idx_20] = min(abs( T_plasma_20 - NarrowTimeFrameTo*day2min ));
    [~,idx_50] = min(abs( T_plasma_50 - NarrowTimeFrameTo*day2min ));
    figure
    yyaxis left
    plot(  T_plasma_100(1:idx_100) ,  E_drug_100(  C_plasma_100(1:idx_100)) ,'LineStyle' , '--','LineWidth',2)
    hold on
    plot(  T_plasma_20(1:idx_20) ,  E_drug_20(  C_plasma_20(1:idx_20)) ,'LineStyle' , '-', 'LineWidth',2)
    hold on 
    plot(  T_plasma_50(1:idx_50) ,  E_drug_50(  C_plasma_50(1:idx_50)) ,'LineStyle' , '-.', 'LineWidth',2)

    ylabel("effects(units)")
    hold on
    yyaxis right
    plot( T_plasma_100(1:idx_100) , C_plasma_100(1:idx_100), 'LineStyle' , '--',  'LineWidth',2    )
    hold on
    plot( T_plasma_20(1:idx_20) , C_plasma_20(1:idx_20), 'LineStyle' , '-', 'LineWidth',2    )
    hold on
    plot( T_plasma_50(1:idx_50) , C_plasma_50(1:idx_50), 'LineStyle' , '-.', 'LineWidth',2    )

    hold off
    ylabel("Plasma concentration(\mug/mL) ")
    hold off
    xticks( 0: 60:  NarrowTimeFrameTo*day2min  )
    xticklabels(  num2cell( 0 :1:NarrowTimeFrameTo*24 )  )
    xlabel("Time(hours)")
    legend(["drug time-effect 100 mg/kg","drug time-effect 20 mg/kg","drug time-effect 50 mg/kg",...
        "plasma concentration 100 mg/kg","plasma concentration 20 mg/kg","plasma concentration 50 mg/kg"],'FontSize',10,'Location','best')
    title("5FU time-effect and time-concentration curve(20,50,100 mg/kg)")
    grid on


%The reason why 100mg/kg with schedule in line with data exhibits over
%inhibitory effects is that the AUC is larger than the one revealed in
%data.

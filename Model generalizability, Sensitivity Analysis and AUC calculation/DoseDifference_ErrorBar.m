clear
close all
addpath("D:\pythonworkshop\matlab\PKPD_5FU\compartmental analysis\DataSet_TGI")
TableProp_NormalizeTV_flag =1;
%The arrangement of the models follow the orders :
 %1. 100 mg/kg,weekly for 4 times  (Colon 38)
 %2. 50 mg/kg, daily for 4 times (Colon 26)
 %3. 20 mg/kg , twice a week for 4 weeks (HT 29)
 %4. 20 mg/kg , three times a week for 4 weeks (HT 29)
 %5. 5 mg/kg, daily for 9 times (SW620)  .
 %6. 20 mg/kg , twice a week for 4 weeks (SW620)
 %7. 20 mg/kg , every the other day for 2 weeks (SW620)
 %8. 50mg/kg ,every the other day ,total 6 times (LS174T sensitive) 
 %9. 50 mg/kg, every the other day , 6 times (LS174T resistant)
 %order of 20mg/kg twice a week CellLine =  ["SW620","HT29","SW480\beta 6"];
 
 delete  ErrorBarTable.xls 
 ErrorBarTable = "ErrorBarTable.xls"; 
 Error_bar_col  = cell(9,1);
 Data_20mgkg_TwiceAWeek_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'20mgkg_TwiceAweek.xls' );
 Data_20mgkg_TwiceAWeek_HT29 = cell2mat( Data_20mgkg_TwiceAWeek_cell(2));
 Error_bar_col {3} =  [  Data_20mgkg_TwiceAWeek_HT29(:,4) Data_20mgkg_TwiceAWeek_HT29(:,5)] ;
 Data_20mgkg_TwiceAWeek_SW620 = cell2mat( Data_20mgkg_TwiceAWeek_cell(1));
 Error_bar_col {6} = [  Data_20mgkg_TwiceAWeek_SW620(:,4) Data_20mgkg_TwiceAWeek_SW620(:,5)];
 %5 daily
 Treatement_error = [0.2150    0.2087    0.3665    0.4910    0.3751    0.3313    0.7134    1.0437    0.62];
 control_error = [0.3896    0.2980    0.3335    0.4938    0.4131    0.8720    1.2900    2.5168    2.7934];
 Error_bar_col {5} = [Treatement_error'  control_error' ];
 % 20 mg/kg , three times a week for 4 weeks (HT 29)
 Data_20mgkg_ThreeAWeek_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'20mgkg_3timesAweek_3weeks.xls');
 Data_20mgkg_ThreeAWeek_HT29 = cell2mat( Data_20mgkg_ThreeAWeek_cell(1));
 Error_bar_col {4} = [   Data_20mgkg_ThreeAWeek_HT29(:,4)   Data_20mgkg_ThreeAWeek_HT29(:,5)];%//
 
 % 20 every the other day SW620
 Data_20mgkg_EveryTheOtherDay_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'20mgkg_EveryTheOtherDay.xls');
 Data_20mgkg_EveryTheOtherDay_SW620 = cell2mat( Data_20mgkg_EveryTheOtherDay_cell(1));
 Error_bar_col {7} = [  Data_20mgkg_EveryTheOtherDay_SW620(:,4)-Data_20mgkg_EveryTheOtherDay_SW620(:,2)   Data_20mgkg_EveryTheOtherDay_SW620(:,5)-Data_20mgkg_EveryTheOtherDay_SW620(:,3)];
 % 50 daily
 Data_50_daily_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'50mgkg_daily_percentage.xls');
 Data_50_daily = cell2mat(Data_50_daily_cell(1));
 Error_bar_col {2}= [ Data_50_daily(:,4)-Data_50_daily(:,2)  Data_50_daily(:,5)-Data_50_daily(:,3)];
 % 50  sensitive
 Data_50_sensitive_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'LS174T_sensitive_5FU_50EveryTheOtherDay.xls');
 Data_50_sensitive = cell2mat(Data_50_sensitive_cell(1));
 Error_bar_col {8} = [Data_50_sensitive(:,4)-Data_50_sensitive(:,2)  Data_50_sensitive(:,5)-Data_50_sensitive(:,3)];
 % 50  resistant
 Data_50_resistant_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'LS174T_resistant_5FU_50EveryTheOtherDay.xls');
 Data_50_resistant = cell2mat(Data_50_resistant_cell(1));
 Error_bar_col {9} = [Data_50_resistant(:,4)-Data_50_resistant(:,2) Data_50_resistant(:,5)-Data_50_resistant(:,3)];
 % 100 weekly
 Data_100_weekly_cell = DoseDiff_TableProcessing (TableProp_NormalizeTV_flag,'TV_100mgkgweekly.xls');
 Data_100_weekly = cell2mat(Data_100_weekly_cell(1));
 Error_bar_col {1} = [Data_100_weekly(:,4)-Data_100_weekly(:,2) Data_100_weekly(:,5)-Data_100_weekly(:,3)];
 for i = 1:9
     ErrorBar = [Error_bar_col{i}];
     writematrix( ErrorBar,  ErrorBarTable,  'Sheet',num2str(i)  );
 end
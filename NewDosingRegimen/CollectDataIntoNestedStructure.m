function AllModel_info   = CollectDataIntoNestedStructure 
%index into a structure that is nested with another strucure.
%parent structure: model
%Model info  example
%      Model_info = {
%         100, "weekly injection for 4 weeks(Colon 38)", [45*day2min 3  7*day2min 92.9936306*plasma_UnitConversion], [497.88706  10.18805 0.68483 5.21749  0.144/day2min  4.8595  0.0432/day2min], 125*mm3tocm3 ,'TV_100mgkgweekly.mat';
%    };
%---------------------------------------------
%The arrangement of the models follow the orders :
 %100 mg/kg,weekly for 4 times  (Colon 38)
 %50 mg/kg, daily for 4 times (Colon 26)
 %20 mg/kg , twice a week for 3 weeks (HT 29)
 %20 mg/kg , three times a week for 4 weeks (HT 29)
 %5 mg/kg, daily for 9 times (SW620)
 %20 mg/kg , twice a week for 3 weeks (SW620)
% 20 mg/kg , every the other day for 2 weeks (SW620)
 %50mg/kg ,every the other day ,total 6 times (LS174T sensitive) 
 %50 mg/kg, every the other day , 6 times (LS174T resistant)
%---------------------------------------------------- 
UnitConversion = [  10^6/130.077, 24*60,10^-3 ,7*24*60   ];
UnitConversion_cell = num2cell(UnitConversion  );
[plasma_UnitConversion, day2min  , mm3tocm3  , week2min  ] = UnitConversion_cell{:};

%T_lag, IC_50, E_max_damage, EC_50_damage, lambda_g, P_max,lambda_d
load NineTGI_Para.mat  MeanPara_array
AllPara = MeanPara_array;
Seed_Model_DisplayLabel_pool = ["weekly injection for 4 weeks(Colon 38)"
    "daily injection, 4 times(Colon 26)"
    "Twice a week for 3 weeks(HT29)"
    "Three times a week for 3 weeks(HT29)"
    "daily for 9 days(SW620)"
    "Twice a week for 3 weeks(SW620)"
    "Every other day for 2 weeks(SW620)"
    "Every  other day 6 times(sensitive LS174T)"
    "Every other day 6 times(resistant LS174T)"
    ];
TV_initial_col = [125 125 100 100  56.76  0.0565/mm3tocm3  71.3 125 125  ] *mm3tocm3;
Dose_info_col = [45*day2min 3  7*day2min 	 92.9936306
     10080  3     day2min  93/2 
    3*7*day2min	3*2-1	7*day2min/2	 93/5
    3*7*day2min	 3*3-1	7*day2min/3 	 93/5 
    9*day2min  8 day2min 	 92.9936306 /20
    3*7*day2min+5*day2min	3*2-1	7*day2min/2	 93/5 
    2*7*day2min    6     2*day2min  	 93/5 
    15*60*24 5 48*60  93/2 
    15*60*24 5 48*60  93/2 
    
    ];
FeaturedDose_col  = [100, 50, 20,20,5,20,20,50,50];
DataSource_col = {'TV_100mgkgweekly.mat'
    '50mgkg_daily.mat'
    '20mgkg_TwiceAweek_HT29.mat'
    '20mgkg_ThreeTimesAweek_HT29.mat'
    '5mgkg_daily_SW620.mat'
    '20mgkg_TwiceAweek_Data1_SW620.mat'
    '20mgkg_EveryTheOtherDay_SW620.mat'
    '50mgkg_EveryTheOtherDay_LS174T_sensitive.mat'
    '50mgkg_EveryTheOtherDay_LS174T_resistant.mat'
    
    };
UpstreamKinetics_col= {
    'Time_Anabolites_DSB_100_Weekly_4times.mat'
    'Time_Anabolites_DSB_50_Daily_4times.mat'
    'Time_Anabolites_DSB_20_TwiceAweek_3weeks.mat'
    'Time_Anabolites_DSB_20_ThreeTimesAWeek_3weeks.mat'
    'Time_Anabolites_DSB_5_daily_9times.mat'
    'Time_Anabolites_DSB_20_TwiceAweek_3times.mat'
    'Time_Anabolites_DSB_20_EveryTheOtherDay_3weeks.mat'
    'Time_Anabolites_DSB_50_EveryTheOtherDay_6times_sensitive.mat'
    'Time_Anabolites_DSB_50_EveryTheOtherDay_6times_resistant.mat'
    };
AllModel_info =  [];
%concatenate structure arrays using the [] operator
ErrorBarTable = "ErrorBarTable.xls"; 
for i = 1:9
    FeaturedDose = FeaturedDose_col(i);
    Label =  Seed_Model_DisplayLabel_pool (i);
    Dose_regime  = Dose_info_col(i,:);
    Para = AllPara(i,:);
    InitialTv = TV_initial_col(i);
    DataSource = fullfile('DataSet_TGI' , DataSource_col{i});
    UpstreamKinetics = fullfile( 'TimeCourses_DSB_Anabolites'  ,UpstreamKinetics_col{i});
    ErrorBar = readtable( ErrorBarTable , 'Sheet',num2str(i)   );
    ErrorBar_array = table2array(ErrorBar); 
    Output_struct  = DataCollection(FeaturedDose, Label,Dose_regime, Para,InitialTv,DataSource,UpstreamKinetics, ErrorBar_array);
    AllModel_info  =  [AllModel_info  Output_struct];
end
 
    function Output_struct  = DataCollection(FeaturedDose, Label,Dose_regime, Para,InitialTv,DataSource,UpstreamKinetics, ErrorBar_array)
    %Create a nonscalar structure that contains several fields.
    Data = importdata(DataSource);
    Dose_regime(1)  = max(Dose_regime(1) ,Data(end,1) ); 
    field1 = 'FeaturedDose';  value1 = FeaturedDose;
    field2 = 'Label';  value2 = Label;
    field3 = 'Dose_regime';  value3 = Dose_regime ;
    field4 = 'Para';  value4 = Para;
    field5 = 'InitialTv';  value5 = InitialTv;
    field6 = 'DataSource';  value6 = DataSource;
    field7 = 'UpstreamKinetics';  value7 = UpstreamKinetics;
    field8 = 'ErrorBar';  value8 =  ErrorBar_array;

    s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8 );
    Output_struct = s;
    end
end
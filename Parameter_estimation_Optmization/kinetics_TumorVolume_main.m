clear  
close all
addpath('../DataSet_TGI')
Dose = 50; %change
estimation_flag =1; %change 
mm3tocm3 = 10^-3; %10^-3
plasma_UnitConversion = 10^6/130.077;  %nug/mL to pmol/mL
day2min = 24*60;

if Dose == 5
    theta_fitted_Tumor_ControlGroupPara = [ ];
elseif Dose == 20
    theta_fitted_Tumor_ControlGroupPara = [ ];
elseif Dose == 50
    theta_fitted_Tumor_ControlGroupPara = [];
elseif Dose == 100
    theta_fitted_Tumor_ControlGroupPara = [ 0.144/day2min  4.8595  0.0432/day2min];  %day^-1   cm^3  day^-1
end
if estimation_flag== 0
    theta0 = [] ;
    theta_fitted_Tumor_known = [
         
% 580.589661618
% 0.000282501
% 1.883650695
% 546.052116600
% 0.000296805
% 3.009205428
% 441.095755798
% 0.000281016
% 1.874742382
 0.000178194  
 1.189253116
 0.000025640];
else   %parameters need to be estimated
%"T_{tumor}","IC_{50}","k_{death}", 
    theta_fitted_Tumor_known  =[];
    theta0 = [
    0.000270951
1.926389378
0.000058033
 ] ;
end
CurveFitting_info = [100  0.75 1.25 2];
%theta_fitted_Tumor_known_ExtraForPlot = []; %vertical
theta_fitted_Tumor = TumorVolume_estimation(estimation_flag, theta_fitted_Tumor_known ,theta0 , theta_fitted_Tumor_ControlGroupPara,...
    CurveFitting_info,  Dose );

function  theta_fitted_Tumor = TumorVolume_estimation(estimation_flag ,theta_fitted_Tumor_known ,theta0 , theta_fitted_Tumor_ControlGroupPara,...
     CurveFitting_info, Dose )
day2min = 24*60;
mm3tocm3 = 10^-3; %10^-3
Weektomin = 7*day2min;
plasma_UnitConversion = 10^6/130.077; 

if Dose ==100 
    InitialTV_Treatment =125*mm3tocm3; % mm^3
%     Data_array = importdata('TV_100mgkgweekly.mat');
%     Data_array(:,1) =  Data_array(:,1)*day2min;
    Data_array = importdata('100mgkg_resistance_Colon26.mat');
    Dose_info =  [45*day2min 3  7*day2min 	 92.9936306*plasma_UnitConversion];
    Dose_DisplayName= '100mg/kg, weekly injection for 4 weeks'  ;
elseif Dose == 50
    InitialTV_Treatment  = 125*mm3tocm3; 
%     Data_array = importdata('50mgkg_daily.mat');
%     Dose_info = [    Data_array(end,1)	  3     day2min  93/2*plasma_UnitConversion]; %weekly
%     Dose_DisplayName= '50mg/kg, daily 4 times' ;
    Dose_info = [15*60*24 5 48*60  93/2*plasma_UnitConversion];
    Data_array = importdata('50mgkg_EveryTheOtherDay_LS174T_resistant.mat');
    Dose_DisplayName= '50mg/kg,every the other day, cell line: LS174T(resistant)' ;
elseif Dose ==20
    InitialTV_Treatment  =71.3*mm3tocm3 ;    %    71.3 for every the other day  ;100 For twice a day
%%20mgkg_TwiceAweek
%     Data_array = importdata('20mgkg_TwiceAweek_Data1_SW620.mat');%one sheet, one cell 
%     Dose_info = [    3*7*day2min+5*day2min	3*2-1	7*day2min/2	 93/5*plasma_UnitConversion]; 
%     Dose_DisplayName= '20mg/kg, injection twice a week for 3 weeks(SW620)' ;
%%20mgkg_ThreeTimesAweek
%    Data_array = importdata('20mgkg_ThreeTimesAweek_HT29.mat');%one sheet, one cell 
%    Dose_info = [   3*7*day2min	 3*3-1	7*day2min/3 	 93/5*plasma_UnitConversion]; 
%    Dose_DisplayName= '20mg/kg, injection three times a week for 3 weeks,cell line:HT29' ;
   Data_array = importdata('20mgkg_EveryTheOtherDay_SW620.mat');%one sheet, one cell 
   Dose_info = [    2*7*day2min    6     2*day2min  	 93/5*plasma_UnitConversion]; 
   Dose_DisplayName= '20mg/kg, injection every the other day for 2 weeks,cell line:SW620' ;
   
elseif Dose == 5 
   InitialTV_Treatment  = 56.76*mm3tocm3;  
   %InitialTV_NoTreatment  = 46.85*mm3tocm3;  
   Data_array = importdata('5mgkg_daily.mat');%one sheet, one cell 
   Dose_info = [       9*day2min  8 day2min 	 92.9936306*plasma_UnitConversion/20]; 
   Dose_DisplayName= '5mg/kg, injection daily for 9 days,cell line: SW620' ;
end
Dose_info_cell = num2cell( Dose_info);
[Duration, Interval_TotalNum,DoseFrequency , ~ ] = Dose_info_cell{:};

t_data_TV = Data_array(:,1);
NonNaN_Rowindex =  find(all(~isnan(Data_array), 2) ) ; 
Relative_data_Treatment = Data_array(:,2);
t_data_TV_NoTreatment =  Data_array(NonNaN_Rowindex,1);
Relative_data_NoTreatment = Data_array(NonNaN_Rowindex,3);

para_num_Tumor = size(theta0,1);
CurveFitting_info_cell = num2cell(CurveFitting_info);
[population_size,  LB,  UB,  MaxIter ] =  CurveFitting_info_cell{:};
Iter_num = 1; %initialize
theta_collection = zeros(MaxIter+1,para_num_Tumor+1);
RSS_collection = zeros(MaxIter,1);
Optimizeoption_flag = 3; %fmincon
multistart_Startpoint_flag =1; %use CustomizeStartPoint to genertate starting points1
myfun_Tumor = @kinetics_TumorVolume_estimation;
myfun_Tumor_plot = @kinetics_TvPlot_DoseDifference;
parameter_label_Control = ["\lambda_{g}","p_{max}","\lambda_{d}"];
parameter_label = ["T_{lag}","IC_{50}", "E_{max,damage}","EC_{50,damage}", parameter_label_Control ];
parameter_unit_label = [" min", " pmol/mg" ," units"," units", "min^{-1}","cm^3", " min^{-1}"]; %total 6 parameters
num_para = size(parameter_unit_label,2);
elongatedTime = 28*day2min; %37*day2min;
 
Dose_info_all = [
  %45*day2min	3	7*day2min	 93/20*plasma_UnitConversion
    %45*day2min	3	7*day2min	 93/10*plasma_UnitConversion
    45*day2min	3	7*day2min	 93/5*plasma_UnitConversion
    %45*day2min	3	7*day2min	 93/2*plasma_UnitConversion
   % 45*day2min	3	7*day2min	 93*plasma_UnitConversion
    45*day2min	3	7*day2min	 0
    ];%Duration, Interval_TotalNum,DoseFrequency ,Dose
%Dose_info_all (:,1:3) = ones(size(Dose_info_all,1),1) * Dose_info(1:3);
Dose_info_all = [ Dose_info ;   45*day2min	3	7*day2min	 0];
Dose_info_all (:,1)  = ones(size(Dose_info_all,1),1)* (Duration + elongatedTime);
Dose_info_all_display = append("Model: ", ["5 mg/kg","10 mg/kg","20 mg/kg","50 mg/kg"]);

if estimation_flag ==1
    if multistart_Startpoint_flag == 0 %RandomStart
        theta0 = rand(para_num_Tumor,1);
        theta_collection(1,1:1:para_num_Tumor) = theta0;
        disp('initial theta');
        disp(theta0);
        theta_collection(Iter_num,1:1:para_num_Tumor ) = theta0;
        while Iter_num <= MaxIter
            start_point_flag = 0;%random start
            [theta_g,RSS] = kinetics5_para_estimation(t_tumor,c_Tumor, theta0, myfun_Tumor, para_num_Tumor, Optimizeoption_flag, start_point_flag,[]) ;
            theta0 = theta_g;
            theta_collection(Iter_num+1,1:1:para_num_Tumor) = theta_g;
            theta_collection(Iter_num+1,para_num_Tumor+1) = RSS;
            RSS_collection(Iter_num) = RSS;
            fprintf('\t\tRSS in %dth iteration\t\t = %8.5f\n', Iter_num, RSS);
            Iter_num = Iter_num+1;
            
        end
        
    end
    
    if multistart_Startpoint_flag == 1 % CustomizeStart
        start_point_flag = 1;
        %fix the parameters for growth kinetics and
        %set the larger initial values of IC50 larger than others
        theta0 = theta0'; %theta0 is a horizonal vector
        init_para = theta0;
        disp('initial theta');
        disp(theta0);
        theta_collection(Iter_num,1:1:para_num_Tumor ) = theta0;
        while Iter_num <= MaxIter
            %use the best fitted theta of the previous iteration as the initial value of parameters in the current iteration
            para_initGuess_modified = initGuess_generater(init_para,population_size,LB,UB, para_num_Tumor,[],[]);
            start_para = para_initGuess_modified;
            [theta_g,RSS] = kinetics5_para_estimation([],[], init_para,myfun_Tumor, para_num_Tumor, Optimizeoption_flag, start_point_flag,start_para) ;
            init_para = theta_g; %1*para_num matrx
            theta_collection(Iter_num+1,1:1:para_num_Tumor ) = theta_g;
            theta_collection(Iter_num+1,para_num_Tumor +1) = RSS;
            RSS_collection(Iter_num) = RSS;
            for k1 = 1:length( theta_g)
                fprintf(1, '\t\tTheta(%d) = %8.9f\n', k1,theta_g(k1))
            end
            fprintf('\t\tRSS in %dth iteration\t = %8.9f\n', Iter_num, RSS);
            [c, ceq] = AUCe_NonlinearConstraint(theta_g);
            fprintf('\t\t Nonlinear constraints\t = %8.9f\n', c(1)  );
            pause(2);
            Iter_num = Iter_num+1;
            %theta_g_plot = [theta_g(1:4)  theta_fitted_Tumor_ControlGroupPara(1)   theta_g(5)   theta_fitted_Tumor_ControlGroupPara(2)   ]; %theta_g is horizontal
         % 20 three times a week
               theta_g_plot = [893.53059    0.00445 0.24860 0.26015  theta_g]   ; 
              % 5m/kg daily
           % theta_g_plot = [893.53059 0.00445 0.24860 0.26015  theta_g ]
%            theta_g_plot = theta_g;
            % 20 every the other day
%            theta_g_plot = [  893.53059  0.00445 0.24860 0.26015 theta_g ];
           % 50 every the other day
%             theta_g_plot = [893.53059 theta_g  0.000249874 7.717727070 0.000004283];
            [Tfit_Absolute_NoTreatment, Cfit_Absolute_NoTreatment]  =  myfun_Tumor_plot(theta_g_plot, InitialTV_Treatment , Dose_info_all(end,:));

            %plot
            figure 
            plot(t_data_TV,  Relative_data_Treatment ,'o', 'MarkerSize',7,'DisplayName', 'data');%data
            hold on
            plot(Tfit_Absolute_NoTreatment, Cfit_Absolute_NoTreatment(:,1)/InitialTV_Treatment , 'LineWidth',2);%model
            hold on
            plot(t_data_TV_NoTreatment,Relative_data_NoTreatment, '^', 'MarkerSize',8, 'DisplayName','Control group TV data '); %data
            hold on
            for i = 1:size(Dose_info_all,1)-1 %no no treatment
            	[Tfit,Cfit] =  myfun_Tumor_plot(  theta_g_plot,InitialTV_Treatment , Dose_info_all(i,:)) ;    
                plot(Tfit,  Cfit(:,1)/InitialTV_Treatment , 'LineWidth',2, 'DisplayName',Dose_info_all_display(i));%model 
            end
            hold off
            grid
            xticks( 0:Weektomin: fix( Dose_info_all(1,1) /Weektomin)*Weektomin )
            xticklabels(  num2cell( 0 :fix( Dose_info_all(1,1)  /Weektomin )  ) )
            xlabel('Time(weeks)')
            ylabel('Relative tumor volume')
            title("Tumor growth" + newline+ strjoin(append( parameter_label,"=", string(theta_g_plot(1: num_para )         ))    ))
            legend('Location','southeast')
            drawnow;    
            
            figure
            [Tfit,Cfit] =  myfun_Tumor_plot(  theta_g_plot,InitialTV_Treatment , Dose_info_all(end-1,:)) ; 
            plot(Tfit,  Cfit(:,2), 'LineWidth',2);%model
            grid
            xticks( 0: DoseFrequency: fix( Dose_info_all(1,1)/DoseFrequency )*DoseFrequency )
            xticklabels(  num2cell( 0 :fix( Dose_info_all(1,1)/DoseFrequency )  ) )
            xlabel('Time(weeks)')
            ylabel('tumor volume')
            title("Damamged cells")
            legend('Location','southeast')
            drawnow;
            
        end
    end
    theta_fitted_Tumor =  theta_collection ( find( RSS_collection == min( RSS_collection ))+1, 1:1:para_num_Tumor);
    RSS_fitted_Tumor = min(RSS_collection);
    for k1 = 1:length( theta_fitted_Tumor)
        fprintf(1, '\t\tTheta(%d) = %8.9f\n', k1,theta_fitted_Tumor(k1))
    end
    fprintf('\t\tRSS\t\t = %8.9f\n',   RSS_fitted_Tumor)
    
elseif estimation_flag == 0
    theta_fitted_Tumor = theta_fitted_Tumor_known ;
    theta_fitted_Tumor = transpose(theta_fitted_Tumor(:)); %column to row vector 
end
%time curve
 
 theta_fitted_Tumor_plot  =  [ 893.53059  0.00445 0.24860 0.26015 theta_fitted_Tumor ]; 
  [Tfit_Absolute_NoTreatment, Cfit_Absolute_NoTreatment]  =  myfun_Tumor_plot(theta_fitted_Tumor_plot,InitialTV_Treatment, Dose_info_all(end,:));

%plot           
figure

if Dose == 5
    Treatement_error = [0.2150    0.2087    0.3665    0.4910    0.3751    0.3313    0.7134    1.0437    0.62];
    control_error = [0.3896    0.2980    0.3335    0.4938    0.4131    0.8720    1.2900    2.5168    2.7934];
    errorbar(t_data_TV, Relative_data_Treatment,Treatement_error,'o', 'MarkerSize',7,'MarkerFaceColor',[0 0.4470 0.7410], 'DisplayName', 'Treated group TV data');%data
    hold on
    errorbar(t_data_TV_NoTreatment,Relative_data_NoTreatment, control_error,'^', 'MarkerSize',7, 'DisplayName','Control group TV data '); %data
    hold on
else
    plot(t_data_TV, Relative_data_Treatment,'o', 'MarkerSize',7,'MarkerFaceColor',[0 0.4470 0.7410], 'DisplayName', 'Treated group TV data');%data
    hold on
    plot(t_data_TV_NoTreatment,Relative_data_NoTreatment, '^', 'MarkerSize',7, 'DisplayName','Control group TV data '); %data
    hold on
end
plot(Tfit_Absolute_NoTreatment, Cfit_Absolute_NoTreatment(:,1)/InitialTV_Treatment, 'LineWidth',2, 'DisplayName',"Simulated untreated TV with "...
    +newline+strjoin( append ( parameter_label_Control ,"=",  string(theta_fitted_Tumor_plot(end-2:end)) ,parameter_unit_label(end-2:end)    )));%model
hold on
for i = 1:size(Dose_info_all,1)-1 %no no treatment
    [Tfit,Cfit] =  myfun_Tumor_plot(theta_fitted_Tumor_plot, InitialTV_Treatment, Dose_info_all(i,:)) ;
    plot(Tfit,  Cfit(:,1)/InitialTV_Treatment, 'LineWidth',2, 'DisplayName',Dose_info_all_display(i));%model
end
hold on
plot( 0: DoseFrequency : Interval_TotalNum*DoseFrequency  , zeros( Interval_TotalNum+1,1   ),   '*' ,'MarkerSize',8, 'DisplayName' ,  Dose_DisplayName )
hold off
grid
xticks( 0:Weektomin: fix( Dose_info_all(1,1) /Weektomin)*Weektomin )
xticklabels(  num2cell( 0 :fix( Dose_info_all(1,1)  /Weektomin )  ) )
xlabel('Time(weeks)')
ylabel(' Relative tumor volume (fold change) ')
title("Normalized Tumor volume based on the starting value" ...
    +"up to" + num2str(fix( (Duration+ elongatedTime)/ 7/day2min  ) ) + 'weeks'+ ...
    newline +"model fit to :" +  Dose_DisplayName +...
    newline+ strjoin( append( parameter_label,"=", string(theta_fitted_Tumor_plot(1: num_para) ), parameter_unit_label,",")     ) , 'FontSize' , 13   )
legend('FontSize' , 12,'Location','southeast')



figure
[Tfit,Cfit] =  myfun_Tumor_plot( theta_fitted_Tumor_plot,InitialTV_Treatment , Dose_info_all(end-1,:)) ;
plot(Tfit,  Cfit(:,2), 'LineWidth',2);%model
hold on
plot( 0: DoseFrequency : Interval_TotalNum*DoseFrequency  , zeros( Interval_TotalNum+1,1   ),   '*' ,'MarkerSize',8, 'DisplayName' ,  Dose_DisplayName )
hold off
grid
xticks( 0:Weektomin: fix( Dose_info_all(1,1) /Weektomin)*Weektomin )
xticklabels(  num2cell( 0 :fix( Dose_info_all(1,1)  /Weektomin )  ) )
xlabel('Time(weeks)')
ylabel('tumor volume(cm^3)')
title("Damamged cell population")
legend('Location','southeast')

%caculate SSR
NonNaN_Rowindex = find( all(~isnan(Data_array), 2) ); 
Cfit_Treatment_q = interp1(Tfit,   Cfit(:,1) ,t_data_TV , 'PCHIP');
Cfit_control_q   =  interp1(Tfit_Absolute_NoTreatment ,   Cfit_Absolute_NoTreatment(:,1)  ,t_data_TV(NonNaN_Rowindex) , 'PCHIP');   
Ysim = [ Cfit_Treatment_q/InitialTV_Treatment    - Relative_data_Treatment ;...
    Cfit_control_q/InitialTV_Treatment  - Relative_data_NoTreatment(NonNaN_Rowindex)];
SSR =  sum( Ysim.^2);
disp(SSR);


Dose_info(1) = Dose_info(1)+    elongatedTime; 
MarkerEdgeColor = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.4940 0.1840 0.5560]];
swap =0;
if Dose == 50  
    %Plot the simulation result for sensitive and resistance cell lines in
    %the same figure
    theta_fitted_LS174T  = [      893.53059   50.051548422 0.835170629 0.040987103   0.000269765  5.826570136 0.000019946
    893.53059  69.380333945 0.526513177 49.468042148 0.000249874 7.717727070 0.000004283 ]; 
    %The control group are all the sensitive one.
    if swap  ==1
      theta_fitted_LS174T(:,end-2:end) =repmat ( theta_fitted_LS174T(1,end-2:end),2,1   );
    end
    Data_array_sensitive = importdata('50mgkg_EveryTheOtherDay_LS174T_sensitive.mat');
    Data_array_resistant = importdata('50mgkg_EveryTheOtherDay_LS174T_resistant.mat');
    Dose_DisplayName =  [ "50mg/kg,every the other day, cell line: LS174T(sensitive)" ,"50mg/kg,every the other day, cell line: LS174T(resistant) "];
    figure 
    plot(Data_array_sensitive(:,1),  Data_array_resistant(:,2),'o', 'MarkerSize',7,'MarkerEdgeColor',MarkerEdgeColor(2,:), 'DisplayName', 'LS174T resistant cell line treated group data');%data
    hold on
    plot(Data_array_sensitive(:,1) , Data_array_resistant(:,3), '^', 'MarkerSize',7,'MarkerEdgeColor',MarkerEdgeColor(2,:), 'DisplayName','LS174T resistant cell line control group data'); %data
    hold on
    plot(Data_array_sensitive(:,1),  Data_array_sensitive(:,2),'o', 'MarkerSize',7,'MarkerEdgeColor',MarkerEdgeColor(1,:), 'DisplayName', 'LS174T sensitive cell line treated group data');%data
    hold on
    plot(Data_array_sensitive(:,1) , Data_array_sensitive(:,3), '^', 'MarkerSize',7,'MarkerEdgeColor',MarkerEdgeColor(1,:), 'DisplayName','LS174T sensitive cell line control group data'); %data
    hold on
    plot( 0: DoseFrequency : Interval_TotalNum*DoseFrequency  , zeros( Interval_TotalNum+1,1   ),   '*' ,'MarkerSize',8, 'DisplayName' , "injection every the other day" )
    hold on 
    for i = 1:2 %plot treated group
        [Tfit,Cfit] =  myfun_Tumor_plot(theta_fitted_LS174T(i,:), InitialTV_Treatment, Dose_info) ;
        plot(Tfit,  Cfit(:,1)/InitialTV_Treatment, 'LineWidth',2, 'DisplayName', Dose_DisplayName(i) ,'Color' ,MarkerEdgeColor(i,:));%model
    end
    hold on
    if swap == 0
        for i = 1:2 %plot control group
            [Tfit,Cfit] =  myfun_Tumor_plot(theta_fitted_LS174T(i,:), InitialTV_Treatment,Dose_info_all(2,:) ) ;
            plot(Tfit,  Cfit(:,1)/InitialTV_Treatment, 'LineWidth',2, 'Color',MarkerEdgeColor(i,:) , 'DisplayName', strcat(Dose_DisplayName(i)," control group"));%model
        end
        title("Normalized Tumor volume based on the starting value" ...
            +"up to" + num2str(fix( (Duration+ elongatedTime)/ 7/day2min  ) ) + 'weeks'+ ...
            newline+"Model for sensitive LS174T:"  +strjoin( append( parameter_label,"=", string(theta_fitted_LS174T(1,1: num_para) ), parameter_unit_label,",")     ) +newline+...
            "Model for resistant LS174T:" +strjoin( append( parameter_label,"=", string(theta_fitted_LS174T(2,1: num_para) ), parameter_unit_label,",")     ), ...
            'FontSize' , 11.5 )
    else
        [Tfit,Cfit] =  myfun_Tumor_plot(theta_fitted_LS174T(1,:), InitialTV_Treatment,Dose_info_all(2,:) ) ;
        plot(Tfit,  Cfit(:,1)/InitialTV_Treatment, 'LineWidth',2, 'Color',MarkerEdgeColor(1,:) , 'DisplayName', strcat(Dose_DisplayName(1)," control group"));%model
        
        title("Normalized Tumor volume based on the starting value" ...
            +"up to" + num2str(fix( (Duration+ elongatedTime)/ 7/day2min  ) ) + 'weeks'+ ...
            newline+ "Control group: sensitive LS174T"...
            + newline +"Model for sensitive LS174T:"  +strjoin( append( parameter_label,"=", string(theta_fitted_LS174T(1,1: num_para) ), parameter_unit_label,",")     ) +newline+...
            "Model for resistant LS174T:" +strjoin( append( parameter_label,"=", string(theta_fitted_LS174T(2,1: num_para) ), parameter_unit_label,",")     ), ...
            'FontSize' , 11.5 )
    end
    grid
    xticks( 0:Weektomin: fix( Dose_info_all(1,1) /Weektomin)*Weektomin )
    xticklabels(  num2cell( 0 :fix( Dose_info_all(1,1)  /Weektomin )  ) )
    xlabel('Time(weeks)')
    ylabel(' Relative tumor volume (fold change) ')
    legend('FontSize' , 12,'Location','best')
    
    



end
end











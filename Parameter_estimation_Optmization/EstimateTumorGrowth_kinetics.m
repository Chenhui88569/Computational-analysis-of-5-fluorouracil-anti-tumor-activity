function  obj_output   = EstimateTumorGrowth_kinetics(theta)
day2min = 24*60;
mm3tocm3 = 10^-3;
Dose = 5 ;

if Dose ==100 
    Tv_initial  = 125*mm3tocm3; % mm^3 
    Data_array = importdata('TV_100mgkgweekly.mat');
    Data_array(:,1)  = Data_array(:,1)*day2min;  
elseif Dose ==20
    Tv_initial  = 100*mm3tocm3; % 20mg/kg gets its own starting value
    Data_array = importdata('20mgkg_TwiceAweek.mat' ); %absolute
elseif Dose ==50
    Tv_initial  = 125*mm3tocm3; % mm^3  
    Data_array = importdata('50mgkg_daily.mat' ); %relative
elseif Dose ==5
    Tv_initial  = 56.76*mm3tocm3;  
    Data_array = importdata('5mgkg_daily.mat');%one sheet, one cell 
 

end
Data_array( any( isnan(Data_array), 2) ,: ) = []; %remove the row 
t_data_TV = Data_array(:,1);
 
Relative_data_TV = Data_array(:,3);%control group
 
c0_Tv  = Tv_initial ;
%NormalTVGrowth_TimePoint =  t_data_TV(end-2:end) - t_data_TV(end-2); %last four points 
[Tv,Cv]=ode15s(@fun_TV,t_data_TV , c0_Tv ); %time span
Ysim =   Cv./Tv_initial   -  Relative_data_TV;%best 
obj_output = Ysim;
% SSR =  sum( Ysim.^2);
% ObjFunction = @(num_datapoints, SSR, num_para) num_datapoints*log(SSR/num_datapoints)+2*(num_para+1)*num_datapoints...
%     /(num_datapoints-num_para-2);
% num_datapoints = size(Ysim,1);
% num_para = 3;
% obj_output =   ObjFunction(num_datapoints,SSR,num_para);


    function diff = fun_TV(t,c)
        %DSB_q = interp1(PlasmaDSB_TimePoint, DSB_TimeCourse , t, 'PCHIP'); 
        dcdt = zeros(1,1);
        lambda_g = theta(1); 
        P_max = theta(2);
        lambda_d  = theta(3);
        dcdt(1)= lambda_g*c(1)*(1- c(1)/P_max) - lambda_d*c(1);        
        diff = dcdt;
    end

end


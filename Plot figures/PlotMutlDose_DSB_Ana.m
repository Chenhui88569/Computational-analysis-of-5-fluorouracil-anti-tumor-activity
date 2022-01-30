clear 
close all

UnitConversion = 10^6/130.077;  %nug/mL to pmol/mL
day2min = 24*60;
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
Dose_info(:,end) = [ 1/20; 1/10;1/10;1/5 ;1/5; 1/2;1/5;1 ;0 ]*Dose_info(end-1,end);
Dose_info_display = append( "Model:", [
    "5 mg/kg, daily for 9 days" ,...
    "10 mg/kg, 3 times a week for 3 weeks" , "10 mg/kg, twice a week for 5 weeks", ...
    "20 mg/kg, every the other day for 2 weeks" , "20 mg/kg, twice a week for 3 weeks",...
    "50 mg/kg, daily, 4 times",...
    "20 mg/kg, 3 times a week for 3 weeks",...
    "100mg/kg, Weekly schedule for 4 weeks" ," control group" ]) ;
for i = 1:size(Dose_info,1)
    Plot_Intermediate(Dose_info(i,:),Dose_info_display(i) );
end

function Plot_Intermediate(Dose_info,sup_title)
%Duration, Interval_TotalNum, DoseFrequency ,Dose
Dose_info_cell = num2cell(Dose_info);
Weektomin = 7*24*60;
[Duration, Interval_TotalNum,DoseFrequency ,Dose ] = Dose_info_cell{:};
COL = kinetics_uptoDSB_alternative(Dose_info);
%title_col =  [ "Plasma" ,"Peripheral", "Interstitial Fluid", "Cellular" ,"Metabolites", "FRNA", "FDNA", "dUMP",...
%  "TS-FdUMP complex","dNTP", "DSB" ,"%free TS"];
title_col =  [ "Metabolites","DSB"];
%Legend_Newposition = [0.010073968719506,0.814332939592635,0.091536459724108,0.042470526056402];
if Dose == 0 % control case
    f = figure;
    subplot(1,2, 1      )
    p = plot(COL(:,1), COL(:,2), 'LineWidth',2  );
    title(  title_col(1)  )
    xticks( 0:Weektomin: fix( Duration/Weektomin)*Weektomin )
    xticklabels(  num2cell( 0 :fix( Duration/Weektomin )  ) )
    grid
    
    
    subplot(1,2, 2      )
    p= plot(COL(:,1), COL(:,3), 'LineWidth',2  );
    xticks( 0:Weektomin: fix( Duration/Weektomin)*Weektomin )
    xticklabels(  num2cell( 0 :fix( Duration/Weektomin )  ) )
    grid
    xlabel('Time(weeks)')
    title(  title_col(2)  )
    %programmatically create a overall legend
    %place the legend of the last subplot to the new position.
    % legend( p(end),"Model", 'Position',Legend_Newposition,'FontSize',11.5 );
    sgtitle( sup_title)
else
    f = figure;
    subplot(1,2, 1      )
    p = plot(COL(:,1), COL(:,2), 'LineWidth',2  );
    hold on
    p2 = plot( 0: DoseFrequency:  Interval_TotalNum*DoseFrequency , zeros( Interval_TotalNum+1,1   ),   '*' );
    hold off
    title(  title_col(1)  )
    xticks( 0:Weektomin: fix( Duration/Weektomin)*Weektomin )
    xticklabels(  num2cell( 0 :fix( Duration/Weektomin )  ) )
    grid
    
    subplot(1,2, 2      )
    p = plot(COL(:,1), COL(:,3), 'LineWidth',2  );
    hold on
    p2 = plot( 0: DoseFrequency:  Interval_TotalNum*DoseFrequency , zeros( Interval_TotalNum+1,1   ),   '*' );
    hold off
    title(  title_col(2)  )
    xticks( 0:Weektomin: fix( Duration/Weektomin)*Weektomin )
    xticklabels(  num2cell( 0 :fix( Duration/Weektomin )  ) )
    grid
    
    
    %programmatically create a overall legend 
    %place the legend of the last subplots to the new position.
 %   legend( [p1(end),p2(end)],["Model","Repeated injection"],  'Position',Legend_Newposition,'FontSize',11  );
    sgtitle( sup_title)
    
end
newtitle = strrep(sup_title,',','_');
newtitle = strrep(newtitle,' ','');
newtitle = strrep(newtitle,':','_');
newtitle = strrep(newtitle,'/','');
saveas(f,  fullfile('../plot/Intermediates_DSB_ana_MultipleDoses' ,strcat(newtitle ,'.png'))   );

end
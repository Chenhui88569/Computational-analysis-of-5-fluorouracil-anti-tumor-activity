function c_dNTP_measure_pertubation = dNTP_data_processing(Current_dNTP_M)
Table_dir = 'DataSet_all.xlsx';
[T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir) ;
c_dNTP_array = cell2mat(  C_allsheets_output_cell(end-2 )  );
t_dNTP_array = cell2mat(  T_allsheets_output_cell(end-2)  );
c_dNTP_M_cell = cell(size(c_dNTP_array, 2),1); 
for i = 1:size(c_dNTP_array , 2) %column
     c_temp  = c_dNTP_array(1:end-1, i)/100; %convert percentage to fold change 
     c_temp(:,2) = c_temp(:,1);
     dNTP_pertubation_measures = sqrt((log(c_temp)).^2);
    c_dNTP_M_cell{i} = dNTP_pertubation_measures;
end
%cell2mat
%order:dGTP	dCTP	dATP	dTTP
c_dGTP_M = cell2mat(c_dNTP_M_cell(1));
c_dCTP_M = cell2mat(c_dNTP_M_cell(2));
c_dATP_M = cell2mat(c_dNTP_M_cell(3));
c_dTTP_M = cell2mat(c_dNTP_M_cell(4));
c_dNTP_control = c_dNTP_array(end,:)*0.2;
c_dNTP_control_total  = sum(c_dNTP_control ,'all');  %Control :(pmoles/mug DNA)
c_dNTP_percent_control = c_dNTP_array(1:size(t_dNTP_array,1),:);
c_dNTP_absolute = c_dNTP_percent_control*c_dNTP_control'/100 ; %multiply according to dimension
c_dNTP_total_percent_control  = c_dNTP_absolute/c_dNTP_control_total;
c_dNTP_total_pertubation_measures =  sqrt((log(c_dNTP_total_percent_control)).^2);
c_dNTP_total_pertubation_measures(:,2) = c_dNTP_total_pertubation_measures(:,1);
if strcmp(Current_dNTP_M,'alldNTP')
    c_dNTP_measure_pertubation = c_dNTP_total_pertubation_measures;
end

if strcmp(Current_dNTP_M,'dATP')
    c_dNTP_measure_pertubation =c_dATP_M;
    
end

if strcmp(Current_dNTP_M,'dCTP')
    c_dNTP_measure_pertubation =c_dCTP_M;
end

if strcmp(Current_dNTP_M,'dGTP')
    c_dNTP_measure_pertubation =c_dGTP_M;
end


if strcmp(Current_dNTP_M,'dTTP')
    c_dNTP_measure_pertubation =c_dTTP_M;
    
if strcmp(Current_dNTP_M,'All') 
    c_dNTP_measure_pertubation = cell(size(c_dNTP_array, 2)+1,1);
    c_dNTP_measure_pertubation{1} = c_dATP_M;
    c_dNTP_measure_pertubation{2} =c_dCTP_M;
    c_dNTP_measure_pertubation{3} =c_dGTP_M;
    c_dNTP_measure_pertubation{4} =c_dTTP_M;
    c_dNTP_measure_pertubation{5} =c_dNTP_total_pertubation_measures;
    
    
    
    
end
end


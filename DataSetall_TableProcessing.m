function  [T_allsheets_output_cell , C_allsheets_output_cell] = DataSetall_TableProcessing(Table_dir)
%one dataset, one sheet
Data_SheetName = sheetnames(Table_dir);
T_allsheets_output_cell = cell(size(Data_SheetName ));
C_allsheets_output_cell = cell(size(Data_SheetName ));
h2min = 60;
day2min = 24*60;
for i  = 1:size(Data_SheetName,1)
    Data_table = readtable( Table_dir, 'Sheet',Data_SheetName(i)   );
    Data_array = table2array(Data_table);
    Data_array_output = Data_array; 
    if contains(Data_SheetName(i),'(h)'  )
         Data_array_output(:,1)  =  Data_array_output(:,1) *h2min;
    elseif contains(Data_SheetName(i),'(day)'  )
         Data_array_output(:,1)  =  Data_array_output(:,1) *day2min;
    end
    t_temp = Data_array_output(:,1);
    T_allsheets_output_cell{i} =  t_temp(~isnan(t_temp)) ;
    C_allsheets_output_cell{i} = Data_array_output(:,2:end);
end 
end


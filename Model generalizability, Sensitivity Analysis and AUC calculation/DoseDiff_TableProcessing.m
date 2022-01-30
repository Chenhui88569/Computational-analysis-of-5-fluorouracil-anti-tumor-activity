function  Data_allsheets_output_cell = DoseDiff_TableProcessing(NormalizeTV_flag,Table_dir)

Data_SheetName = sheetnames(Table_dir);
Data_allsheets_output_cell =cell(size(Data_SheetName ));
week2min = 7*60*24;
day2min  = 60*24;
for i  = 1:size(Data_SheetName,1)
    Data_table = readtable( Table_dir, 'Sheet',Data_SheetName(i)   );
    Data_array = table2array(Data_table);
    Data_array_output = Data_array; 
    if  Data_array_output(1,1) ~= 0 %make sure that the data points start at 0 day
        Data_array_output(:,1) = Data_array(:,1) - Data_array(1,1); 
    end
    %There are various time scales of the Tumor volume data from different efficacy
    %studies. Some of them are measured regularly and the time points are labeled as day.
    %The others are specifed as weeks and are saved as seperate tables or seperate sheets within one excel file. 
    %The names of those excel file contain 'LabelWeek'.
    [uniq, ~, ~] = unique(round(Data_array_output(:,1),0 ) ); %to distinguish the intergral time points from the decimal time points
    if  contains( Table_dir,"LabelWeek")  %format the timescale to weeks 
        if  size(uniq) == size(Data_array_output(:,1)) %Time points scaled as weeks can be viewed as a series of integers
            Data_array_output(:,1) =  round(Data_array_output(:,1),0 );
        else %Time points scaled as weeks are decimal.  
            Data_array_output(:,1) =  round(Data_array_output(:,1),2 );%Therefore, they need to be rounded to 2 digits to the right of the decimal point
        end
        Data_array_output(:,1) = Data_array_output(:,1)*week2min;
    else % the timescale is already in days
        if  size(uniq) == size(Data_array_output(:,1)) 
            Data_array_output(:,1) =  round(Data_array_output(:,1),0 );%time points are supposed to be a series of integers
        else 
            Data_array_output(:,1) =  round(Data_array_output(:,1),2 );
        end
        Data_array_output(:,1) = Data_array_output(:,1)*day2min ;
    end
    %column 2th is the tumor volume data of treated groups. Column 3th is tumor
    %volume of untreated groups. 
    if  NormalizeTV_flag ==1
        Data_array_output(:,[2,3]) = Data_array(:,[2,3])./Data_array(1,[2,3]); %normalized TV data of intial value.
        if size(Data_array,2) > 3
            Data_array_output(:,4) = Data_array(:,4)./Data_array(1,2); %normalized treated error bar
            Data_array_output(:,5) = Data_array(:,5)./Data_array(1,3); %normalized control error bar 
        end
        Data_allsheets_output_cell{i} =  Data_array_output ;
        
    elseif NormalizeTV_flag ==0%don't need to be normalized
        Data_allsheets_output_cell{i} =  Data_array_output ;
    end
end

%I once considered outer joining all the PD data using the Time as key
%variable. Non-matched key variable in both joined tables are marked as
%'NaN'. But the result return didn't meet my expectation. The
%shoud-be-matched key variable didn't show up in the same row. So this
%idea was dropped. 

%Data_table_output = array2table(Data_table);
% Data_table_allVar  = 1:width(  Data_table);
% Data_table_allVar_minus1 = 1:width(  Data_table)-1;
% Data_table_NewVarname  = append('Data', string( Data_table_allVar_minus1  ));
% Data_table_NewVarname = ["Time"  Data_table_NewVarname ];
% Data_table_output = renamevars( Data_table, Data_table_allVar  ,Data_table_NewVarname );


end


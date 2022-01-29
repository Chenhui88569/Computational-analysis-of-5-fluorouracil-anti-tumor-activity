function output_legend  = CompareWithOriginalValue(Concerned_Para_variation, Concerned_Para_original ,ParaLabel) 
%Concerned_Para_variation and ParaLabel should be vertical vector, The dimension of  Concerned_Para_original doesn't matter  
%ParaLabel_output and string_para_value_ouput should have the same dimension after the for loops ended, so that the two vectors can be concatenated by append fucntion.
%notice the difference between string array and cell array of character vectors
%strjoin(C) linking the elements of C with a space between consecutive elements. no matter C is vertical or horizontal

%The save space in the legend, the untested parameters are eliminated.
Num_para = size(  Concerned_Para_variation,1  );
string_para_value  = cell(Num_para,1  );
for  n = 1:Num_para
    temp1 = Concerned_Para_variation(n);
    temp = regexprep(char(vpa(  temp1)),'([0-9]+\.[0-9]+)','${num2str(str2num($1),''%.3e'')}');
    string_para_value{n} =    temp ;
end

% Concerned_Para_variation_round = round(Concerned_Para_variation,3);
% string_para_value =  cellstr(num2str(Concerned_Para_variation_round)); %to convert matrix (double) to cell array of string(cell arrays of character vectors)
Num_para = size(  Concerned_Para_variation,1  );
ParaLabel = ParaLabel'; %the dimension of ParaLabel should be the same as Concerned_Para_variation
ParaLabel_output = cell( size( ParaLabel  ));  %defined as cell array of characters string( size( ParaLabel    ));
string_para_value_ouput =  cell( size(  string_para_value));%string( size( string_para_value  )) ; 
for i = 1: Num_para
    if Concerned_Para_variation(i) == Concerned_Para_original(i) %original
        ParaLabel_output{i} = 'NaN';
      % ParaLabel_output{i} = ParaLabel{i}  ;
        string_para_value_ouput{i} = 'NaN'; 
       %string_para_value_ouput{i} = append( string_para_value{i}, '(-)'); 
    elseif Concerned_Para_variation(i) > Concerned_Para_original(i)
       FoldChange = round(Concerned_Para_variation(i ) /Concerned_Para_original(i),1) ; 
       string_para_value_ouput{i} = append( string_para_value{i},'(',num2str( FoldChange ),' fold $\uparrow$)' ) ;
       ParaLabel_output{i} = ParaLabel{i}  ;
    elseif Concerned_Para_variation(i) < Concerned_Para_original(i)
       FoldChange = round(Concerned_Para_original(i)/Concerned_Para_variation(i),1);
       string_para_value_ouput{i} = append( string_para_value{i},'(',num2str(FoldChange),' fold $\downarrow$)' );
       ParaLabel_output{i} = ParaLabel{i} ;
    end
end

string_para_value_ouput = string_para_value_ouput(~strcmp( string_para_value_ouput,'NaN'    )    );
ParaLabel_output = ParaLabel_output(~strcmp(  ParaLabel_output,'NaN'    ) ); %find non NaN element  
%output_legend = append( ParaLabel_output,'=',string_para_value_ouput,',') ; %cell array
output_legend  = string_para_value_ouput;
%isnan is used to detect nan value in the numerical array eg.A = [nan 1 2
%], isnan(A)
%isnan can't operate on cell array and string array
%ismissing is used to find the missing value(nan) in the string array.

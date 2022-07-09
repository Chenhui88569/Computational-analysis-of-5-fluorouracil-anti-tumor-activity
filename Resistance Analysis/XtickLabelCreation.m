function [t_xtick, t_xlabel] = XtickLabelCreation(t_list, t_list_label)
%the unit of the elements in t_list is min
t_xtick = reshape(t_list,[1,size(t_list,1)]);
t_xlabel = reshape(t_list_label ,[1,size(t_list,1)]) ;
t_xlabel = round(t_xlabel,2);
t_xlabel  = num2cell(t_xlabel,1);
end


function  d_col = get_flag_modelsim(outfile, delta,Current_para, traj_index_col, Para_Actual_value,ParaNeedToBeChnagedAtEach_Step_col)


day2min = 24*60;
T = (5:5:45)*day2min;
Model_Index = 1;
num_traj_OneBatch = length(traj_index_col);
d_col = zeros(num_traj_OneBatch,length(T));
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% FigHandle = FigList(1);
% set(0, 'CurrentFigure', FigHandle);
for n = 1:num_traj_OneBatch
        traj_idx = traj_index_col(n);
        f_A  =0 ;
        f_B = 0;
        Para_Actual_value_elem =  Para_Actual_value(:,:,  traj_idx );
        path =ParaNeedToBeChnagedAtEach_Step_col(:,  traj_idx );
        idx = find( path == Current_para);
        flag = 0;
        try
            %isemptywarning('  ');
            [Tsim_A,Csim_A,f_A] =  kinetics_TGI(  Para_Actual_value_elem(:,idx) , Model_Index,T);
            differences_A = diff(f_A);
            second_diff_A = diff( differences_A );
            [Tsim_B,Csim_B,f_B] =  kinetics_TGI(  Para_Actual_value_elem(:,idx+1) , Model_Index,T);
            differences_B = diff(f_B);
            second_diff_B =  diff( differences_B );
%             if ~isempty( lastwarn )
            %if  sum(differences_A <0) >=1 &&    sum( differences_B<0) >= 1
           if  abs(max(differences_A )) > 5 || abs(max(differences_B ))  >5 || abs(max(second_diff_A  )) > 3|| abs(max(second_diff_B))  >3
                flag = 1;
           end
        catch
            flag =1;
        end
        if flag == 1
            d_col(n,:) = NaN(1,length(T));
            fprintf(   outfile, ' %d th parameter set fails\n' ,     traj_idx);
            continue
        else
%             m1 = f_A(end);%end);
%             m2 = f_B(end);%end);
            fprintf( outfile, ' %d th parameter set succeeds\n' ,     traj_idx);
%             plot(Tsim_A,Csim_A)
%             hold on
%             plot(Tsim_B,Csim_B)
%             drawnow
            d_col(n,:) = abs(f_A(:)' - f_B(:)')/delta;
        end
end
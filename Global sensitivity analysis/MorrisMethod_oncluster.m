clear
close all
dbstop if error
%warning('off','all')
load('Para.mat')
rng('shuffle')
num_Morris  = 3;
mu = Para_Info.Para_Value;
Para_Name =  Para_Info.Para_Name;
num_para = length(mu);
%Each trajecotory has P+1 parameter
%Using five Morris levels (p = 5), each parameter value is increased or decreased by 25% of its standard value (Δ = 1/(p − 1) = 1/4).
%generate 1000 random samples.
num_traj = 5000;
ub_col = zeros(num_para, 1);
lb_col = zeros(num_para,1);
%dist_col = zeros(num_traj,1);
for m = 1:num_para
    ub_col(m) = mu(m)*1.5;
    lb_col(m) = mu(m)*1/2;
end
corr_col  = [0;1/5;2/5;3/5;4/5;1];   %[0,1/,2/3,1];
delta = 3/5;

for iter = 1:num_Morris
    outdir = ['outputfile_Morris' num2str(iter) '/'];
    if ~isfolder(outdir)
        mkdir(outdir);
    else 
        delete(outdir)
    end
    dir_log = dir([outdir '*.log']);
    if ~isempty( [dir_log.name])
        delete([outdir '*.log'])
    end

    Path_In_Hypercube_col = zeros(num_para, num_para + 1, num_traj);
    Para_Actual_value = zeros(num_para, num_para + 1, num_traj );
    ParaNeedToBeChnagedAtEach_Step_col = zeros(num_para, num_traj);
    for n  = 1:num_traj %59 samples points on each path
        % temp = randi(num_para , num_para ,1);
        temp = randperm(num_para);
        ParaNeedToBeChnagedAtEach_Step_col(:,n) =  temp';
        temp2 = rand(num_para,1);
        temp3 = ones(num_para,1);
        temp3(  temp2 >= 0.5) = -1;
        temp4 = zeros(num_para ,num_para+1);
        temp5 =  randi(length(corr_col), num_para,1); %wrong
        %sample_path(:,1)  = [1/2*ones(5,1); 1/2*ones(27,1); corr_col( temp5(32+1:32+18) ); corr_col( temp5(32+19:end) )];
        sample_path(:,1)  = corr_col(  temp5);
        for j = 2: num_para+1 %  number of sample points
            actual_index = j-1;
            order = temp(actual_index);
            e_delta = zeros(num_para,1);
            e_delta(order) = delta*temp3(actual_index);
            new_corr =  sample_path(:,j-1) + e_delta; % deviation from the previous sample point.
            if sum(new_corr>1) == 1
                idx = find(new_corr>1);
                new_corr( idx ) =   sample_path(idx,j-1) - delta;
            elseif sum(new_corr<0) == 1
                idx = find(new_corr<0);
                new_corr(   idx ) =   sample_path(   idx,j-1) + delta;
            end
            sample_path(:,j) = new_corr;
        end
        Path_In_Hypercube_col (:,:,n) =      sample_path;
        Para_Actual_value(:,:,n) = repmat(lb_col, 1, num_para +1) +    sample_path.*repmat(ub_col-lb_col, 1, num_para+1);
        %% time saving
        %     centered_sample_path_norm = vecnorm(sample_path,2,1); %l2 norm
        %     dist = sum(centered_sample_path_norm);
        %     dist_col(n) =  dist;
        %%
    end
    % selected r = 25 samples in a way to maximize their “spread” in the input space.
    % clustering the sample path based on the the geomertric distances between
    % the sample paths and the standard paths.
    
    num_traj_used = 30;
    day2min = 24*60;
    T = (5:5:45)*day2min;
    % [idx,C,sumd,D] = kmedoids(  dist_col,num_cluster) ;
    % mediod_idx = find( ismember( dist_col, C));
    % x_col = 1:1:num_traj;
    % trun_Path_Actual_value =  Path_Actual_value(:,:, mediod_idx);
    % trun_para_path_col = ParaNeedToBeChnagedAtEach_Step_col(:,mediod_idx);
    % figure
    % for n = 1:num_cluster
    %     scatter(x_col(idx==n),dist_col(idx==n))
    %     hold on
    %     plot(mediod_idx(n),dist_col(mediod_idx(n)),'co',...
    %         'MarkerSize',7,'LineWidth',1.5)
    % end
    %num_traj_used*50
    d_col_cells = cell(num_para,1);
    
    % final tumor size as the metric
    
    Model_Index = 1;
    %[Tsim,Csim,f_A] =  kinetics_TGI( mu , Model_Index,T);
    %f1 = figure;
    parpool(num_para)
    parfor k =  1: num_para%1:num_para
        task  = getCurrentTask;
        taskid = task.Id;
        filename = [ outdir,'logfile_' num2str(taskid) '.log'];
        outfile = fopen(  filename , 'a');
        fprintf(   outfile, 'Analysis on %d th parameter, %s\n' ,   k,  Para_Name(k));
        d_col = get_flag_modelsim( outfile, delta, k,1: num_traj_used ,  Para_Actual_value,ParaNeedToBeChnagedAtEach_Step_col) ;
        num_non_nan = sum(~isnan(d_col(:,1) )); % the number of two consecutive trajectories getting nan
        
        num_iter  = 1;
        fprintf( outfile, '%d th try, number of success in total %d\n' ,   num_iter,    num_non_nan);
        temp_col = zeros(num_iter,1);
        temp_col(1) = num_traj_used;
        while  num_non_nan < num_traj_used  && num_iter <=  6
            new_num = (num_traj_used- num_non_nan)*7;
            num_iter  =  num_iter +1;
            temp_col( num_iter) = new_num;
            %         if num_iter == 1
            %             d_col_new = get_flag_modelsim(outfile, delta, k,   temp_col(1)+1 : cumsum( temp_col),...
            %             Para_Actual_value,ParaNeedToBeChnagedAtEach_Step_col) ;
            %         else
            % new_set =  num_traj_used+  (num_iter-1)* new_num+1 : num_traj_used+ num_iter * new_num;
            cumsum_new_set = cumsum( temp_col);
            new_set =  cumsum_new_set(end-1)+1 : cumsum_new_set(end);
            d_col_new = get_flag_modelsim(outfile, delta, k, new_set,   Para_Actual_value,ParaNeedToBeChnagedAtEach_Step_col) ;
            %         end
            d_col = [ d_col; d_col_new];
            num_non_nan  = sum(~isnan(d_col(:,1) ));
            
            fprintf(   outfile,'%d th try, number of success %d\n',   num_iter,    num_non_nan);
        end
        if num_non_nan < num_traj_used
            fprintf(   outfile,'fail on %s, number of success in total %d\n', Para_Name(k),    num_non_nan);
        else
            fprintf(   outfile,'Success on %s, number of success in total %d\n', Para_Name(k),    num_non_nan);
        end
        d_col_cells{k} = d_col(~isnan(d_col(:,1)),:) ;
        fclose(outfile);
    end
    delete(gcp('nocreate'))
    save([outdir  'Morris_col_cells' num2str(iter) '.mat'],'d_col_cells');
end
%         flag_B = 0;

%         try
%             waring('');
%             [Tsim,Csim,f_B] =  kinetics_TGI(  Para_Actual_value_elem(:,idx+1) , Model_Index,T);
%             if ~isempty(lastwarn)
%                 flag_B = 1;
%             end
%         catch
%             flag_B =1;
%         end
%         if flag_B == 1
%             d_col(k,n) = nan;
%             continue
%         else
%             m2 = f_B(end);
%         end
%         clear f_A f_B

    %d_col_no_nan = d_col(k,~isnan(d_col (k,:)    ));
    

% d_col_no_nan = d_col(~isnan(d_col));
% d_mean_col = mean(d_col_no_nan , 2);
% d_mean_col_norm = d_mean_col/max(d_mean_col);
% centered_d_col = d_col  - repmat(d_mean_col, 1,num_cluster);
% std = vecnorm(centered_d_col,2,2)/sqrt(num_cluster-1);
% std_norm = std/max(std);
% figure
% h = bar([d_mean_col_norm std_norm ]);
% set(h, {'DisplayName'}, {'\mu^*','\sigma'}')
% save('MorrisResult.mat','d_mean_col_norm','std_norm');
    
    
    
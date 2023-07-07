function  [Para_col,  acc_rate_col,    S_curr]  = MyRAM(  T, num_para, Para_col ,S_curr,myfun,Var_name, outputfile_dir,Model_Index  )
curr_T  = 1 ;

if ~isempty(Model_Index)
    Target_curr = myfun( Para_col(1,:), Model_Index )  ;
else 
    Target_curr = myfun( Para_col(1,:))  ;
end
num_figure = ceil(num_para/6);
for i = 1:  num_figure
    figure;
end
% %Var_name  =  {'V_{max}' , 'Q_{21}' ,    'V_1',  'V_2',  'K_m','\sigma'};
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
[~,idx] = sort( [FigList.Number]);
FigList  = FigList(idx);

alpha_bar = 0.234; %

acc_rate_col = zeros(1, T);
gamma= 0.2;
%outputfile_dir = 'MCMC_cellular/';
if ~isfolder(outputfile_dir)
    mkdir(outputfile_dir)
else 
    rmdir(outputfile_dir,'s');
     mkdir(outputfile_dir)
end
outfile = fopen( [outputfile_dir  'outlog.log'],'a');
while  curr_T  <= T 
    curr_para = Para_col(curr_T,:);
    curr_para = curr_para';
    U = randn(num_para,1);
    next_para = curr_para +  S_curr*U;
    while  ~prod( next_para > 0) 
        U = randn(num_para,1);
        next_para = curr_para +  S_curr*U;
    end
    if ~isempty(Model_Index)
        Target_next = myfun( next_para, Model_Index )  ;
    else
        Target_next = myfun( next_para)  ;
    end
    %Target_next =  myfun(next_para,Model_Index)  ;
    %Target_next  =  Target_next * PriorFun_PKPD(next_para );
    u = rand;
    alpha_i =  min([1,   exp( log(Target_next)  - log(Target_curr))       ]) ;
    if  u<   alpha_i  % Target_next/Target_curr
        Para_col(curr_T+1,:)  =  next_para'; 
        acc_rate_col(curr_T) = 1;
        fprintf(outfile, 'accept , %d \r', curr_T );
        Target_curr =     Target_next ;
    else
        Para_col(curr_T+1,:)  =   curr_para; 
        Target_curr = Target_curr;
    end
    if rem(curr_T,50) == 0 %100*floor(curr_T/100) = curr_T  
        fprintf(outfile,'current iteration:  %d \r ', curr_T);
    end
    
    factor  =    U*U'./vecnorm( U,2)^2;
    temp_mat = S_curr*(  eye(num_para) + min(1,curr_T^(-gamma))*(   alpha_i- alpha_bar  )*factor  )*S_curr';
    S_curr = chol(temp_mat, 'lower');
    flag = 0 ;
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        set(0, 'CurrentFigure', FigHandle);
        for j = 1:6
            subplot(2,3, j)
            idx = (iFig-1)*6+j;
            plot(Para_col(1: curr_T, idx),'k-','LineWidth',1.5)
            xlabel('iteration')
            %axis([1+T*fix(curr_T/T),T+T*fix(curr_T/T),0,15]);
            grid off
            title(Var_name(idx))
            drawnow
            if idx == num_para
                flag  = 1;
                break
            end
        end
        if flag == 1
            break
        end 
    end
    curr_T  = curr_T  +1 ;
end
end
  

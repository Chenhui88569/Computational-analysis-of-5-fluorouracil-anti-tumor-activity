function [theta_g,RSS,manymins] = kinetics5_para_estimation(t,c,theta0,myfun,para_num,OptimizeSolver_flag,multistart_Startpoint_flag,start_para) 
% remember to remove many times
%lb = 1.0e-5*ones(para_num,1);
lb =  [0; 
    zeros(para_num-1,1)];
%lb(3) = 5.0e-6 ; 
ub = inf*ones(para_num,1);
%ub  =   10*ones(para_num,1);

 %In the case of 100mg/kg, the carrying capacity of treated group is larger
 %that the untreatd group. Therefore, the upper limit is posed on \chi_max
x_data = t;
y_data = c;
%FitOptions = optimset('Display', 'iter-detailed', 'Algorithm', 'levenberg-marquardt');
%options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','ScaleProblem','jacobian');
%train plasama without options 
%train tumor with options
%theta0 = start_para(1,:);
MutiStartPoint_num = 200;

if OptimizeSolver_flag ==1 
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','ScaleProblem','jacobian');
   [thetafitted,errorfitted]=lsqcurvefit(myfun,theta0,x_data,y_data,lb,ub,options);
    problem = createOptimProblem('lsqcurvefit','x0',theta0,'objective',myfun,...
    'lb',lb,'ub',ub,'xdata',x_data,'ydata',y_data);
    ms = MultiStart('UseParallel',true); 
    if multistart_Startpoint_flag ==0 
        tpoints = RandomStartPointSet('NumStartPoints',MutiStartPoint_num);
        [theta_fit,RSS_fit] = run(ms,problem,tpoints); %theta_fit: the best point found, fval contains the squared norm of the residual.
        theta_g =  theta_fit;
        RSS = RSS_fit;
    end
    if multistart_Startpoint_flag ==1
       tpoints = CustomStartPointSet(start_para);
       [theta_fit,RSS_fit] = run(ms,problem,tpoints);
       theta_g =  theta_fit;
       RSS = RSS_fit;
    end
 
end

if OptimizeSolver_flag == 0 %Optimizeoption:MaxIterations
    %options = optimoptions('lsqcurvefit','MaxIterations',100,'Display','iter','MaxFunctionEvaluations',2000);
    options = optimoptions('lsqcurvefit','Display','iter');
    [thetafitted,errorfitted]=lsqcurvefit(myfun,theta0,x_data,y_data,lb,ub,options);
    problem = createOptimProblem('lsqcurvefit','x0',theta0,'objective',myfun,...
    'lb',lb,'ub',ub,'xdata',x_data,'ydata',y_data);
    ms = MultiStart('UseParallel',true,'Display','iter'); 
    %gs = GlobalSearch(ms);%Valid solvers are: fmincon.
    
    %ms = MultiStart;
    if multistart_Startpoint_flag ==0
        tpoints = RandomStartPointSet('NumStartPoints',MutiStartPoint_num);
        [theta_fit,RSS_fit] = run(ms,problem,tpoints);
        theta_g =  theta_fit;
        RSS = RSS_fit;
    end
    if multistart_Startpoint_flag ==1
       tpoints = CustomStartPointSet(start_para);
       [theta_fit,RSS_fit,eflag,output,manymins]  = run(ms,problem,tpoints);
        theta_g =  theta_fit;
        RSS = RSS_fit;
    end
end

if OptimizeSolver_flag ==2 %when the time frames of the species in a certain submodel are not consistent
    options = optimoptions('lsqnonlin','Display','iter');
    [thetafitted,errorfitted]=lsqnonlin(myfun,theta0,lb,ub,options);
    problem = createOptimProblem('lsqnonlin','x0',theta0,'objective',myfun,...
    'lb',lb,'ub',ub);
    ms = MultiStart('UseParallel',true,'Display','iter'); 
    %gs = GlobalSearch(ms);%Valid solvers are: fmincon.
    
    %ms = MultiStart;
    if multistart_Startpoint_flag ==0
        tpoints = RandomStartPointSet('NumStartPoints',MutiStartPoint_num);
        [theta_fit,RSS_fit] = run(ms,problem,tpoints);
        theta_g =  theta_fit;
        RSS = RSS_fit;
    end
    if multistart_Startpoint_flag ==1
       tpoints = CustomStartPointSet(start_para);
        [theta_fit,RSS_fit] = run(ms,problem,tpoints);
        theta_g =  theta_fit;
        RSS = RSS_fit;
    end
  
end


if OptimizeSolver_flag ==3 %when the time frames of the species in a certain submodel are not consistent
   opts = optimoptions('fmincon','Display','iter');
   nonlcon = @AUCe_NonlinearConstraint;
   [thetafitted,errorfitted] = fmincon(myfun,theta0,[],[],[],[],lb,ub,nonlcon,opts);
   problem = createOptimProblem('fmincon','x0',theta0,'objective',myfun,...
    'lb',lb,'ub',ub ,'nonlcon', nonlcon );
   ms = MultiStart('UseParallel',true,'Display','iter'); 
    %gs = GlobalSearch(ms);%Valid solvers are: fmincon.
    
    %ms = MultiStart;
    if multistart_Startpoint_flag ==0
        tpoints = RandomStartPointSet('NumStartPoints',MutiStartPoint_num);
        [theta_fit,RSS_fit] = run(ms,problem,tpoints);
        theta_g =  theta_fit;
        RSS = RSS_fit;
    end
    if multistart_Startpoint_flag ==1
       tpoints = CustomStartPointSet(start_para);
        [theta_fit,RSS_fit] = run(ms,problem,tpoints);
        theta_g =  theta_fit;
        RSS = RSS_fit;
    end
  
end
    


 

end
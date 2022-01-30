function SensitivityAnalysis_TGI(AllModel_info)
addpath('../Resistance')
%1.The baseline parameters is chosen as the ones for model 3, since the model
%has more prominant capability of accomodating different clinical protocols
UnitConversion = [  10^6/130.077, 24*60,10^-3 ,7*24*60   ];
UnitConversion_cell = num2cell(UnitConversion  );
[plasma_UnitConversion, day2min  ,mm3tocm3  ,week2min  ] = UnitConversion_cell{:};
%1. Total number of parameter sets that are drawn out from the normal distribution over a certain range.
%2. The more samples are generated, the 90% confidence interval would be much narrower, the mean of the population time course would be much closer to
%the actual one.
ParaPopulation_size = 800; 
ParaLabel_Control = ["\lambda_{g}","p_{max}","\lambda_{d}"];
ParaLabel = ["T_{lag}","IC_{50}", "E_{max,damage}","EC_{50,damage}", ParaLabel_Control]; 
num_model = 9;
Tv_kinetics_fun  = @kinetics_TvPlot_DoseDifference;
num_para = 7;
% The specified variable appears inside a parfor loop within different indexing expressions. 
% Because the indices are inconsistent across the uses of the array created by the parfor loop, 
% MATLAB sends the entire array to each worker, resulting in high data communication overhead
theta_Tv_temp =  extractfield(AllModel_info,'Para')  ; % the same as theta_Tv_col_1 =   [AllModel_info(:).Para] ;
theta_Tv_col  =  reshape(theta_Tv_temp ,num_para,num_model ); 
theta_Tv_col = theta_Tv_col';
Dose_info_temp  = extractfield(AllModel_info,'Dose_regime');
Dose_info_col  = reshape(Dose_info_temp,4,num_model);
Dose_info_col = Dose_info_col';
TV_initial_col  =  extractfield(AllModel_info,'InitialTv');
C_col = cell(num_para, num_model*2 ); 
for i= 1:num_para
        folderpath = fullfile('..','plot', 'SensitivityAnalysis');
        if ~isfolder(folderpath)
            mkdir(  folderpath);
%         else 
%             delete( fullfile(folderpath, '*')  );
        end
        ParaLabel_removed = strrep(ParaLabel(i),'{','');
        ParaLabel_removed = strrep(ParaLabel_removed,'}','');
        ParaLabel_removed = strrep(ParaLabel_removed,',','');
        C_col{i, 1} = strcat('$',ParaLabel(i),'$');
        path1 = fullfile(folderpath, strcat(ParaLabel_removed,".png" ) ); 
        path2 = fullfile(folderpath, strcat(ParaLabel_removed,".fig" ) ); 
        A = figure('Position', get(0, 'Screensize'));
        for m = 1:num_model
            subplot(3,3,m)
            fprintf("para : %d,Model: %d \n",i,m)
            theta_Tv= theta_Tv_col(m,:); 
            BaseLine_para = theta_Tv; %horizontal vector
            Dose_info = Dose_info_col(m,:);
            elongation = 3*week2min;
            Dose_info(1) = Dose_info(1)+elongation;
            TV_initial =TV_initial_col(m); 
            TumorVolume_TimeSpan  = 0 : 50: Dose_info(1);
            TumorVolume_Timelength =  size(TumorVolume_TimeSpan,2  ); 
            
            TumorVolume_PopulationKinetics =  zeros(ParaPopulation_size, TumorVolume_Timelength  );
            ParaPopulation = repmat(BaseLine_para,ParaPopulation_size,1 ); %baseline population
            UpperBound_fold = 5; % 15% increase = 1.15 fold increase ; 200%increase = 3 fold increase = triple increase.
            LowerBound_fold = 5;
            ub = UpperBound_fold* BaseLine_para(i);
            lb =  1/LowerBound_fold * BaseLine_para(i);
            Population_target =  lb +   ( ub- lb ) .*rand( ParaPopulation_size ,1    ); %vertical vector
            perturbed_idx = randperm(ParaPopulation_size);%shuffle the samples
            ParaPopulation(:,i) = Population_target(perturbed_idx);
            parfor  j =1:ParaPopulation_size
                [Tv,Cfit]  = Tv_kinetics_fun(ParaPopulation(j,:)  ,TV_initial, Dose_info ); %20 mg/kg
                Cfit_Tv = transpose(Cfit(:,1)/TV_initial);
                TumorVolume_PopulationKinetics(j,:) = interp1(Tv ,Cfit_Tv   ,TumorVolume_TimeSpan,  'PCHIP');
            end
            [CI, C_mean, C_MAX, C_MIN, para_MAX, para_MIN ] = Resistance_InferentialAnalysis(...
                TumorVolume_PopulationKinetics, ParaPopulation_size, ParaPopulation    );%CI: confidence interval
            [T_baseline,Cv_baseline] = Tv_kinetics_fun(BaseLine_para ,TV_initial, Dose_info );
            f  = plot( TumorVolume_TimeSpan, [ C_MAX , C_MIN],'-.', TumorVolume_TimeSpan,C_mean ,'--', 'LineWidth',2    );
            hold on
            f1= plot(T_baseline,Cv_baseline(:,1)/TV_initial, 'LineWidth',2    );
            hold on
            [ha hb hc] = shadedplot( TumorVolume_TimeSpan', transpose(CI(:,1)),transpose(CI(:,2)), 'c','none');
            %ha is a 1*2 area array ,  corresponds to two descriptive labels
            hold off
            max_legend   = CompareWithOriginalValue( para_MAX, BaseLine_para,ParaLabel ); %return a cell array of character vector
            min_legend  = CompareWithOriginalValue( para_MIN, BaseLine_para,ParaLabel );
%             legend( [f'], { ['max ', newline ,strjoin(max_legend) ], ...
%                 ['min ',newline,  strjoin( min_legend) ], 'mean'...
%                 },   'Location','southeast')
            grid on
            xlabel('Time(week)')
            ylabel('fold change')
            xticks( 0:week2min: fix( Dose_info /week2min)*week2min )
            xticklabels( (  num2cell( 0:1: fix( Dose_info /week2min) )) )
            title( strcat(" model " , num2str(m) )    )  
            C_col{i, 2*(m-1) + 2} =cell2mat( max_legend);
            C_col{i, 2*(m-1) + 3}= cell2mat( min_legend);
        end
        sgtitle( strcat("Sensitivity analysis on parameter " , ParaLabel(i)  ))
        saveas(gcf, path1)
        savefig(gcf, path2)
end
T = cell2table(C_col)
save Sen_T.mat T
writetable(T,'Sen_T.xls')
end
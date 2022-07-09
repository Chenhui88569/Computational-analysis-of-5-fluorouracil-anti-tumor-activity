# Computational-analysis-of-5-fluorouracil-anti-tumor-activity
The entire model is composed of PK, cellular, and tumor growth inhibition (TGI) sub-models (shown in figures below), the last two of which constitute the PD model. The PK model is used to compute the concentration-time profile, and the TGI model focuses on the kinetics of tumor growth post-treatment.

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gyw1z94b3uj30gl05fgln.jpg" alt="PK" style="zoom:33%;" />
<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gyw1z8ealyj30zo0u0acg.jpg" alt="CellularTumor" width= 700 />


## Resistance analysis
The resistance analysis is discussed in detail in the manucript. The main functions in the folder are for resistance analysis are @Resistance_main_oncluster.m and Resistance_main_plot.m, which run on the cluster using Matlab Parallel computing to save extensive computational time. @Resistance_main_oncluster.m is served to exexcute hte 
  
## Plot figures
Functions for plotting the figures are put in the paper. The scripts are @PlotAllIntermediates.m, @PlotAllTGIModels.m, and @PlotMutIDose_DSB_Ana.m
## parameter estimation
Functions for estimating parameters for each submodel are included. The main functions are kinetics_dNTP_main.m kinetics_DSB_main.m kinetics_PKand3to10_main.m kinetics_TumorVolume_main.m. We conduct the sequential parameter estimation because it is computationally infeasible to estimate all the parameters once and all.
## Model generalizability, Sensitivity Analysis, and AUC calculation
The main function is @DoseDifference main.m used to carry out the three tasks.
## TimeCourses_DSB_Anabolites
The folder includes the precomputed time courses of DSB and 5-FU anabolites using predefined cellular model for each of nine TGI models to speed up the process of estimating parameters in TGI models.

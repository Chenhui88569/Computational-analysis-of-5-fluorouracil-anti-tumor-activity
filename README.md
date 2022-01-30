# Computational-analysis-of-5-fluorouracil-anti-tumor-activity
The entire model is composed of PK, cellular, and tumor growth inhibition (TGI) sub-models (shown in figures below), the last two of which constitute the PD model. The PK model is used to compute the concentration-time profile, and the TGI model focuses on the kinetics of tumor growth post-treatment.





<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gyw1z94b3uj30gl05fgln.jpg" alt="PK" style="zoom:33%;" />

<img src="https://tva1.sinaimg.cn/large/008i3skNgy1gyw1z8ealyj30zo0u0acg.jpg" alt="CellularTumor" style="zoom:33%;" />

### parameter estimation
@EstimateTumorGrowth_kinetics : 
  Input arguement: parameter vector
  function used for estimating parameters in TGI model for the untreated group.    
@EstimateTumorGrowth_plot : 
  Input arguement: parameter vector
  function used for estimating parameters in TGI model for the untreated group.  

# Computational-analysis-of-5-fluorouracil-anti-tumor-activity
The entire model is composed of PK, cellular, and tumor growth inhibition (TGI) sub-models (shown in figures below), the last two of which constitute the PD model. The PK model is used to compute the concentration-time profile, and the TGI model focuses on the kinetics of tumor growth post-treatment.





![PK](https://tva1.sinaimg.cn/large/008i3skNgy1gyw1vi8qiij30gl05fgln.jpg)

![CellularTumor](https://tva1.sinaimg.cn/large/008i3skNgy1gyw1vh03zqj30zo0u0acg.jpg)

### parameter estimation
@EstimateTumorGrowth_kinetics : 
  Input arguement: parameter vector
  function used for estimating parameters in TGI model for the untreated group.    
@EstimateTumorGrowth_plot : 
  Input arguement: parameter vector
  function used for estimating parameters in TGI model for the untreated group.  

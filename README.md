# expgumbel - Two parameter exponentiated Gumbel distribution: properties and estimation with flood data example

__Maintainer__: Enayetur Raheem

Dey, S., Raheem, E., Mukherjee, S. and Ng, H.K.T., 2017. [Two parameter exponentiated Gumbel distribution: properties and estimation with flood data example](https://www.tandfonline.com/doi/abs/10.1080/09720510.2016.1228261). Journal of Statistics and Management Systems, 20(2), pp.197-233.

## Abstract

This article addresses various properties and different methods of estimation of the parameters of exponentiated Gumbel distribution from the frequentist point of view. 

Various mathematical and statistical properties of the exponentiated Gumbel distribution, such as quantiles, moments, conditional moments, hazard rate function, mean residual lifetime, mean deviation about mean and median, entropies and order statistics, are derived. 

We briefly describe different frequentist approaches, namely, maximum likelihood estimation, method of moments, percentile-based estimation method, least squares estimation, method of maximum product of spacings, method of Cramér-von-Mises and methods based on Anderson-Darling statistic. 

Monte Carlo simulations are performed to compare the performance of the estimation methods for small and large samples. 

The application of the model is studied using a flood data example. Bootstrap method was used to obtain bias and standard error of the estimates as well as the percentile confidence intervals. 

Further, confidence regions for the parameters are obtained using likelihood ratio based method. Finally, the concept of return period is used to predict the occurrence of flood in the future.

The paper can be freely downloaded from [Researchgate](https://www.researchgate.net/profile/Sanku-Dey/publication/313965705_Two_Parameter_Exponentiated_Gumbel_Distribution_Properties_and_Estimation_with_Flood_Data/links/58fb79b6aca2723d79d841bc/Two-Parameter-Exponentiated-Gumbel-Distribution-Properties-and-Estimation-with-Flood-Data.pdf)

## What's in this Repository

- `expGubmel_core_functions.R`: [Core functions for parameter estimation in Exponentiated Gubmel distribution](https://github.com/raheems/expgumbel/blob/main/expGumbel_core_functions.R)

  - Random number generation from EG distribution (includes density, distribution, quantilie, hazard functions)
  
  - Frechet Distribution
  - Generalized Extremevalue Distribution
  - Maximum Likelihood Estimation (MLE)
  - Least squares estimation (LSE)
  - Weighted Least Squared Estimation (WLSE)
  - Percentile Estimation (PCE)
  - Method of Maximum Product Spacing (MPS)
  - Cramer von Mises Estimator (CVM)
  - Anderson Darling Estimator (AD)
  - Right-tailed Anderson Darling Estimator
  - Method of Moment Estimator (MME)
  - Issues with tied cases
  
- `expGubmel_simulation.R`: [Monte Carlo Simulation Experiment](https://github.com/raheems/expgumbel/blob/main/expGumbel_simulation.R)

  - includes core functions to simulation
  - summarizing the simulation results


## Reference

Dey, S., Raheem, E., Mukherjee, S. and Ng, H.K.T., 2017. [Two parameter exponentiated Gumbel distribution: properties and estimation with flood data example](https://www.tandfonline.com/doi/abs/10.1080/09720510.2016.1228261). Journal of Statistics and Management Systems, 20(2), pp.197-233.


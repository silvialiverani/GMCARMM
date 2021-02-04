# GMCARMM
Repository for the code to run a GMCARMM (generalised multivariate conditional autoregressive model with multiple membership) for misaligned areal data.

The method is described in the paper available on [arxiv](https://arxiv.org/abs/2004.05334): 

Marco Gramatica, Peter Congdon and Silvia Liverani (2020) *Bayesian modelling for spatially misaligned health areal data: a multiple membership approach*. Journal of the Royal Statistical Society Series C (to appear). Available at [https://arxiv.org/abs/2004.05334](https://arxiv.org/abs/2004.05334) 

The application in the paper concerns diabetes. 

**Abstract**

Diabetes prevalence is on the rise in the UK, and for public health strategy, estimation of relative disease risk and subsequent mapping is important. We consider an application to London data on diabetes prevalence and mortality. In order to improve the estimation of relative risks we analyse jointly prevalence and mortality data to ensure borrowing strength over the two outcomes. The available data involves two spatial frameworks, areas (middle level super output areas, MSOAs), and general practices (GPs) recruiting patients from several areas. This raises a spatial misalignment issue that we deal with by employing the multiple membership principle. Specifically we translate area spatial effects to explain GP practice prevalence according to proportions of GP populations resident in different areas. A sparse implementation in Stan of both the MCAR and GMCAR allows the comparison of these bivariate priors as well as exploring the different implications for the mapping patterns for both outcomes. The necessary causal precedence of diabetes prevalence over mortality allows a specific conditionality assumption in the GMCAR, not always present in the context of disease mapping. 



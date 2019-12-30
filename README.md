# AsymmetricRisk.jl

Traditional measures of risk, such as the variance of a time series or the covariance between two time series, use all available data points (aka information) when computing the statistical measure. However, it is becoming common to consider separately the notions of "upside" and "downside" risks. For example, a large literature has documented that people fear losses more than they like gains of similar magnitude. This package contains implementations of several asymmetric risk measures that have been published in financial and economics research.

_This package is currently under development. Posting to the Julia General Registry will take place once all features listed below are ready._

Univariate measures:
 - [ ] Variance-based (semi-variance)
 - [ ] Lower partial moments (Harlow 1991)
 - [ ] VaR-based (Expected Shortfall / Expected Longrise)
 
Bivariate measures:
 - [ ] Correlation-based (Ang & Chen 2001 + statistical test of Hong, Tu, & Zhou 2007 + one-sided as in Schreindorfer 2019)
 - [ ] Beta-based (Ang, Chen, & Xing 2006)
 - [ ] Entropy-based (Jiang, Wu, & Zhou 2018)
 - [ ] Coskewness (Harvey & Siddique 2000)
 - [ ] Quantile dependence
  



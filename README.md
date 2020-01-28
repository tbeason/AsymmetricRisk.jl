# AsymmetricRisk.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tbeason.github.io/AsymmetricRisk.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://tbeason.github.io/AsymmetricRisk.jl/dev)

[![Build Status](https://travis-ci.org/tbeason/AsymmetricRisk.jl.svg?branch=master)](https://travis-ci.org/tbeason/AsymmetricRisk.jl)
[![codecov.io](http://codecov.io/github/tbeason/AsymmetricRisk.jl/coverage.svg?branch=master)](http://codecov.io/github/tbeason/AsymmetricRisk.jl?branch=master)

Traditional measures of risk, such as the variance of a time series or the covariance between two time series, use all available data points (aka information) when computing the statistical measure. However, it is becoming common to consider separately the notions of "upside" and "downside" risks. For example, a large literature has documented that people fear losses more than they like gains of similar magnitude. This package contains implementations of several asymmetric risk measures that have been published in financial and economics research.

_This package is currently under development._

Univariate measures:
 - [X] Variance-based (semi-variance)
 - [X] Lower partial moments (Harlow 1991)
 - [ ] VaR-based (Expected Shortfall / Expected Longrise)
 
Bivariate measures:
 - [ ] Correlation-based (Ang & Chen 2001 + statistical test of Hong, Tu, & Zhou 2007 + one-sided as in Schreindorfer 2019)
 - [X] Beta-based (Ang, Chen, & Xing 2006)
 - [ ] Entropy-based (Jiang, Wu, & Zhou 2018)
 - [X] Coskewness (Harvey & Siddique 2000)
 - [X] Quantile dependence

## Usage
From the Julia REPL type
```
]add https://github.com/tbeason/AsymmetricRisk.jl.git
```
which will clone the repository.

You can then use the package by typing
```
using AsymmetricRisk
```

See the [docs](https://tbeason.github.io/AsymmetricRisk.jl/stable) for more information.




  



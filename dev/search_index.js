var documenterSearchIndex = {"docs":
[{"location":"univariate/#Univariate-Functions-1","page":"Univariate Functions","title":"Univariate Functions","text":"","category":"section"},{"location":"univariate/#","page":"Univariate Functions","title":"Univariate Functions","text":"Modules = [AsymmetricRisk]\nPages = [\"univariate.jl\"]","category":"page"},{"location":"univariate/#AsymmetricRisk.expectedshortfall-Tuple{AbstractArray{T,1} where T,Float64}","page":"Univariate Functions","title":"AsymmetricRisk.expectedshortfall","text":"expectedshortfall(X::AbstractVector,alpha)\n\nReturns the nonparametric estimate of Expected Shortfall at the alpha-level for X. Assumes X represents payoffs, not losses (which would be -X).\n\n\n\n\n\n","category":"method"},{"location":"univariate/#AsymmetricRisk.expectedshortfall-Tuple{Distributions.Distribution{Distributions.Univariate,Distributions.Continuous},Float64}","page":"Univariate Functions","title":"AsymmetricRisk.expectedshortfall","text":"expectedshortfall(D::ContinuousUnivariateDistribution,alpha)\n\nReturns the Expected Shortfall at the alpha-level for the distribution D. Assumes D is the distribution of payoffs X, not losses (which would be -X).\n\n\n\n\n\n","category":"method"},{"location":"univariate/#AsymmetricRisk.lowerpartialmom-Tuple{Any,Any,Any}","page":"Univariate Functions","title":"AsymmetricRisk.lowerpartialmom","text":"lowerpartialmom(x,t,n)\n\nn-th Lower Partial Moment for vector-like x using target t. \n\nSource: Harlow (1991)\n\n\n\n\n\n","category":"method"},{"location":"univariate/#AsymmetricRisk.semivariance-Tuple{Any}","page":"Univariate Functions","title":"AsymmetricRisk.semivariance","text":"semivariance(x; kind::Symbol=:lower)\n\nReturns the semivariance of vector-like x, which is the variance computed using only the observations below the mean. \n\nDefaults to lower semivariance, but returns upper semivariance if kind=:upper is supplied.\n\n\n\n\n\n","category":"method"},{"location":"univariate/#AsymmetricRisk.upperpartialmom-Tuple{Any,Any,Any}","page":"Univariate Functions","title":"AsymmetricRisk.upperpartialmom","text":"upperpartialmom(x,t,n)\n\nn-th Upper Partial Moment for vector-like x using target t. \n\nSource: Harlow (1991)\n\n\n\n\n\n","category":"method"},{"location":"univariate/#AsymmetricRisk.valueatrisk-Tuple{AbstractArray{T,1} where T,Float64}","page":"Univariate Functions","title":"AsymmetricRisk.valueatrisk","text":"valueatrisk(X::AbstractVector,alpha)\n\nReturns the nonparametric estimate of VaR (value at risk) at the alpha-level for X.\n\n\n\n\n\n","category":"method"},{"location":"univariate/#AsymmetricRisk.valueatrisk-Tuple{Distributions.Distribution{Distributions.Univariate,Distributions.Continuous},Float64}","page":"Univariate Functions","title":"AsymmetricRisk.valueatrisk","text":"valueatrisk(D::ContinuousUnivariateDistribution,alpha)\n\nReturns the VaR (value at risk) at the alpha-level for the distribution D.\n\n\n\n\n\n","category":"method"},{"location":"bivariate/#Bivariate-Functions-1","page":"Bivariate Functions","title":"Bivariate Functions","text":"","category":"section"},{"location":"bivariate/#","page":"Bivariate Functions","title":"Bivariate Functions","text":"Modules = [AsymmetricRisk]\nPages = [\"bivariate.jl\"]","category":"page"},{"location":"bivariate/#AsymmetricRisk.cokurtosis-Tuple{Any,Any}","page":"Bivariate Functions","title":"AsymmetricRisk.cokurtosis","text":"cokurtosis(x,y;kind::Symbol=:asymmetric)\n\nComputes the cokurtosis between vectors x and y. The default is asymmetric cokurtosis.\n\nkind=:asymmetric: fracE(x - Ex)^3 (y-Ey)sigma_x^3 sigma_y\n\nkind=:symmetric: fracE(x - Ex)^2 (y-Ey)^2sigma_x^2 sigma_y^2\n\n\n\n\n\n","category":"method"},{"location":"bivariate/#AsymmetricRisk.coskewness-Tuple{Any,Any}","page":"Bivariate Functions","title":"AsymmetricRisk.coskewness","text":"coskewness(x,y)\n\nComputes the coskewness between vectors x and y. The vector x receives the powers in the expression.\n\nfracE(x - Ex)^2 (y-Ey)sigma_x^2 sigma_y\n\n\n\n\n\n","category":"method"},{"location":"bivariate/#AsymmetricRisk.downsidebeta-Tuple{Any,Any}","page":"Bivariate Functions","title":"AsymmetricRisk.downsidebeta","text":"downsidebeta(x,y;kind=:both)\n\nComputes the downside beta of y on x as defined in Ang & Chen. \n\nkind=:onx allows to condition only on x.\n\n\n\n\n\n","category":"method"},{"location":"bivariate/#AsymmetricRisk.exceedancecor-Tuple{Any,Any,Any}","page":"Bivariate Functions","title":"AsymmetricRisk.exceedancecor","text":"exceedancecor(x,y,c;kind=:both,side=:lower)\n\nComputes exceedance correlations between x and y, using exceedance threshold c. Keyword kind refers to which type of conditioning: :both conditions on x and y simultaneously, while :onx conditions only on x. Keyword side refers to look below (:lower) or above (:upper) the threshold c. Inputs x and y are standardized before conditioning,  so c should be in standard deviation units.\n\n\n\n\n\n","category":"method"},{"location":"bivariate/#AsymmetricRisk.quantiledep-Tuple{Any,Any,Any}","page":"Bivariate Functions","title":"AsymmetricRisk.quantiledep","text":"quantiledep(x,y,q)\n\nQuantile Dependence at the q-th quantile of vector-like data x and y.\n\nLet q_x and q_y be the q-th quantiles of x and y.\n\nReturns Exy  x leq q_x y leq q_y  q if q<=0.5 and (1 -2q + Exy  x leq q_x y leq q_y)  (1-q) otherwise.\n\n\n\n\n\n","category":"method"},{"location":"bivariate/#AsymmetricRisk.upsidebeta-Tuple{Any,Any}","page":"Bivariate Functions","title":"AsymmetricRisk.upsidebeta","text":"upsidebeta(x,y;kind=:both)\n\nComputes the upside beta of y on x as defined in Ang & Chen. \n\nkind=:onx allows to condition only on x.\n\n\n\n\n\n","category":"method"},{"location":"#AsymmetricRisk.jl-Documentation-1","page":"Introduction","title":"AsymmetricRisk.jl Documentation","text":"","category":"section"},{"location":"#","page":"Introduction","title":"Introduction","text":"AsymmetricRisk.jl is a Julia package for computing asymmetric risk measures on univariate and bivariate data.","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"CurrentModule = AsymmetricRisk","category":"page"},{"location":"#","page":"Introduction","title":"Introduction","text":"Pages = [\"univariate.md\",\"bivariate.md\"]\nDepth = 1","category":"page"}]
}

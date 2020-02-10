


@doc raw"""
    semivariance(x; kind::Symbol=:lower)

Returns the semivariance of vector-like `x`, which is the variance
computed using only the observations below the mean. 

Defaults to lower semivariance, but returns upper semivariance if `kind=:upper` is supplied.
"""
function semivariance(x; kind::Symbol=:lower)
    bigN = length(x)
    mu = mean(x)
    if kind==:upper
        tmp = filter(z -> z > mu,x)
    else
        tmp = filter(z -> z <= mu,x)
    end
    return sum((tmp .- mu).^2)/bigN
end








"""
    lowerpartialmom(x,t,n)

`n`-th Lower Partial Moment for vector-like `x` using target `t`. 

Source: Harlow (1991)
"""
function lowerpartialmom(x,t,n)
    bigN = length(x)
    return sum(filter(z->z>=0,t .- x).^n)/bigN
end


"""
    upperpartialmom(x,t,n)

`n`-th Upper Partial Moment for vector-like `x` using target `t`. 

Source: Harlow (1991)
"""
function upperpartialmom(x,t,n)
    bigN = length(x)
    return sum(filter(z->z>=0,x .- t).^n)/bigN
end





"""
    valueatrisk(D::ContinuousUnivariateDistribution,alpha)

Returns the VaR (value at risk) at the `alpha`-level for the distribution `D`.
"""
valueatrisk(D::Distributions.ContinuousUnivariateDistribution,alpha::Float64) = quantile(D,1-alpha)

"""
    valueatrisk(X::AbstractVector,alpha)

Returns the nonparametric estimate of VaR (value at risk) at the `alpha`-level for `X`.
"""
function valueatrisk(X::AbstractVector,alpha::Float64)
    ecdf = EmpiricalCDF()
    append!(ecdf,X)
    sort!(ecdf)
    icdf = finv(ecdf)
    return icdf(1-alpha)
end



"""
    expectedshortfall(D::ContinuousUnivariateDistribution,alpha)

Returns the Expected Shortfall at the `alpha`-level for the distribution `D`.
Assumes `D` is the distribution of payoffs `X`, not losses (which would be `-X`).
"""
function expectedshortfall(D::ContinuousUnivariateDistribution,alpha::Float64)
    return first(quadgk( x -> -(1/alpha)*valueatrisk(D,x),0.0,alpha))
end


"""
    expectedshortfall(X::AbstractVector,alpha)

Returns the nonparametric estimate of Expected Shortfall at the `alpha`-level for `X`.
Assumes `X` represents payoffs, not losses (which would be `-X`).
"""
function expectedshortfall(X::AbstractVector,alpha::Float64)
    ecdf = EmpiricalCDF()
    append!(ecdf,X)
    sort!(ecdf)
    icdf = finv(ecdf)
    return first(quadgk( x -> -(1/alpha)*icdf(1-x),0.0,alpha))
end


# Bertsimas, Lauprete, Samarov (2004 JEDC) good source for ES/VaR/LPM




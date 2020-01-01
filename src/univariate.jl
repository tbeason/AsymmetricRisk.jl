


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





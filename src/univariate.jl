


"""
    lowerpartialmom(x,t,n)

`n`-th Lower Partial Moment for vector-like `x` using target `t`. 

Source: Harlow (1991)
"""
function lowerpartialmom(x,t,n)
    bigN = length(x)
    return sum(filter(y->y>=0,t .- x).^n)/bigN
end





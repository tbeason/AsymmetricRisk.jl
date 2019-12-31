
"""
quantiledep(x,y,q)

Quantile Dependence at the `q`-th quantile of vector-like data `x` and `y`.
"""
function quantiledep(x,y,q)
    qx = quantile(x,q)
    qy = quantile(y,q)
    tmp = mean((x.<=qx).*(y.<=qy))
    if q <= 0.5
        return tmp/q
    else
        return (1-2*q + tmp)/(1-q)
    end
end




#### TODO
# - better API



#### Correlation-based (Ang & Chen 2001 + statistical test of Hong, Tu, & Zhou 2007 + one-sided as in Schreindorfer 2019)
function downsidecor(x,y)
    qs = 0.15:0.025:1
    out = zeros(length(qs))
    for (i,qi) in enumerate(qs)
        qxi = quantile(x,qi)
        loc = x .<= qxi
        tmpxi = x[loc]
        tmpyi = y[loc]
        out[i] = cor(tmpxi,tmpyi)
        #println([length(tmpxi) length(tmpyi) length(tmpxj) length(tmpyj)])
    end
    return out
end

function upsidecor(x,y)
    qs = 0:0.025:0.85
    out = zeros(length(qs))
    for (i,qi) in enumerate(qs)
        qxi = quantile(x,qi)
        loc = x .>= qxi
        tmpxi = x[loc]
        tmpyi = y[loc]
        out[i] = cor(tmpxi,tmpyi)
        #println([length(tmpxi) length(tmpyi) length(tmpxj) length(tmpyj)])
    end
    return out
end

function angchen(x,y)
    qs = 0.15:0.025:0.5
    out = zeros(length(qs)+1)
    for (i,qi) in enumerate(qs)
        qj = 1 - qi
        qxi = quantile(x,qi)
        qyi = quantile(y,qi)
        qxj = quantile(x,qj)
        qyj = quantile(y,qj)
        tmpxi = filter(xx-> xx <=qxi,x)
        tmpyi = filter(yy-> yy <=qyi,y)
        tmpxj = filter(xx-> xx >qxj,x)
        tmpyj = filter(yy-> yy >qyj,y)
        out[i] = cor(tmpxi,tmpyi)
        out[end-(i-1)] = cor(tmpxj,tmpyj)
    end
    return out
end


##### Asymmetry test from Guofu Zhou et al.
standardize(x) = (x .- mean(x)) ./ std(x)
bartlett(z) = abs(z) < 1 ? (1-abs(z)) : zero(eltype(z))
# this only does 1 value of c
function buildEC(x,y,c)
	T = length(x)
	xs = standardize(x)
	ys = standardize(y)
	iL = (xs .< -c) #.& (ys .< -c)
	iU = (xs .> c) #.& (ys .> c)
	lower_xs = xs[iL]
	upper_xs = xs[iU]
	lower_ys = ys[iL]
	upper_ys = ys[iU]
	Tl = length(lower_xs)
	Tu = length(upper_xs)
	lower_cor = cor(lower_xs,lower_ys)
	upper_cor = cor(upper_xs,upper_ys)
	
	# standardize again
	xs[iL] = standardize(lower_xs)
	ys[iL] = standardize(lower_ys)
	xs[iU] = standardize(upper_xs)
	ys[iU] = standardize(upper_ys)
	
	EC = (((T-Tu)/T) .* (xs.*ys .- upper_cor) .* iU) .- (((T-Tl)/T) .* (xs.*ys .- lower_cor) .* iL)
	return (EC=EC,corL=lower_cor,corU=upper_cor)
end

# this wraps it together
# conditioning is standardized x .< -c[i] vs standardized x .> c[i]
# p represents some sort of bandwith? -- lag order?
function asymmetrictest(x,y,c,p)
	T = length(x)
	Nc = length(c)
	
	if Nc == 1
		outs=buildEC(x,y,c)
		ECmat=outs.EC
		uppercors = outs.corU
		lowercors = outs.corL
		S = 0.0
		for L=-p:p
			w = bartlett(L/p)
			g_l = 0.0
			for t=(abs(L)+1):T
				g_l += ECmat[t]*ECmat[t-abs(L)]
			end
			S += w*(1/T)*g_l
		end
	else
		ECmat = zeros(eltype(x),T,Nc)
		uppercors = zeros(eltype(x),Nc)
		lowercors = zeros(eltype(x),Nc)
		for i =1:Nc
			outs=buildEC(x,y,c[i])
			ECmat[:,i]=outs.EC
			uppercors[i] = outs.corU
			lowercors[i] = outs.corL
		end
		S = zeros(eltype(x),Nc,Nc)
		for L=-p:p
			w = bartlett(L/p)
			g_l = zeros(eltype(x),Nc,Nc)
			for ic=1:Nc,jc=1:Nc
				for t=(abs(L)+1):T
					g_l[ic,jc] += ECmat[t,ic]*ECmat[t-abs(L),jc]
				end
			end
			S += w*(1/T).*g_l
		end
	end
	
	J = T * (uppercors - lowercors)' * inv(S) * (uppercors - lowercors)
	pval = pdf(Chisq(Nc),J)
	return (pval=pval,J=J,S=S,corL=lowercors,corU=uppercors,T=T)
end



#### Entropy-based (Jiang, Wu, & Zhou 2018)
quadrantprobabilities(f,c;kwargs...) = (LQP=hcubature(f, [-10;-10], [-c;-c];kwargs...), UQP=hcubature(f, [c;c], [10;10];kwargs...))
function entropy_S(f,c;kwargs...)
    newf(x) = (sqrt(abs(f(x))) - sqrt(abs(f(-x))))^2
    int,se = hcubature(x->newf(x),[c;c],[20;20];kwargs...)
    return (int/2,se)
end
function DOWN_ASY_JWZ(f,c=0;kwargs...)
    lqp,uqp = quadrantprobabilities(f,c;kwargs...)
    S = entropy_S(f,c;kwargs...)
    val = sign(first(lqp)-first(uqp))*first(S)
    return val
end
# FF = kde((standardize(dfnew[!,:GDPg]),standardize(dfnew[!,:NDSg])))
# PF = InterpKDE(FF)
# quadrantprobabilities.(x->pdf(PF,x[1],x[2]),0:0.5:1.5)
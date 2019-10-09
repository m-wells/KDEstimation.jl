@inline function kde(Kh::Kernel{Scaled,D}, X::AbstractVector, ::Direct) where {D<:Distribution}
    return DirectUnivariateKDE(Kh,X)
end

function (kde::DirectUnivariateKDE)(x::Real)
    n = length(kde.X)
    f(t) = kde.Kh(x - t)
    return sum(f, kde.X)/n
end

function (kde::DirectUnivariateKDE)(X::AbstractVector)
    retval = Vector{Float64}(undef,length(X))
    for (i,x) in enumerate(X)
        retval[i] = kde(x)
    end
    return retval
end

function lscv(Kh::Kernel{Scaled,D}, X::AbstractVector, ::Direct) where {D<:Distribution}
    n = length(X)
    retval = 0.0
    for xi in X
        for xj in X
            u = (xi - xj)
            retval += T_H(Kh,u)
        end
    end
    retval /= n
    retval += 2*Kh(0)
    return retval/n
end

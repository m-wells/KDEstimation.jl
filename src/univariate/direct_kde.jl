function kde(Kh::Kernel{Scaled,D}, x::AbstractVector, ::Direct) where {D<:Distribution}
    return DirectUnivariateKDE(Kh,x)
end

function (kde::DirectUnivariateKDE)(x::Real)
    n = length(kde.x)
    f(xᵢ) = kde.Kh(x - xᵢ)
    return sum(xᵢ -> kde.Kh(x - xᵢ), kde.x)/n
end

function (kde::DirectUnivariateKDE)(xs::AbstractVector)
    retval = Vector{Float64}(undef,length(xs))
    for (i,x) in enumerate(xs)
        retval[i] = kde(x)
    end
    return retval
end

function lscv(Kh::Kernel{Scaled,D}, xs::AbstractVector, ::Direct) where {D<:Distribution}
    n = length(xs)
    retval = 0.0
    for xi in xs
        for xj in xs
            u = (xi - xj)
            retval += T_H(Kh,u)
        end
    end
    retval /= n
    retval += 2*Kh(0)
    return retval/n
end

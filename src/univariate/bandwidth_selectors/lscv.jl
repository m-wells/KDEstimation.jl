struct UnivariateLSCV{D,EM,H,T,Res} <: LSCV{D,EM,Res,1}
    eval_method::EM
    h_pre::H            # bandwidths to use for linesearch (before optimize)
    l_pre::Vector{T}    # corresponding lscv values
    res::Res

    function UnivariateLSCV(
        ::Type{D},
        em::EvaluationMethod,
        h_pre::AbstractVector,
        l_pre::AbstractVector,
        res::Optim.OptimizationResults,
    ) where {D<:Distribution}
        return new{
            D,
            typeof(em),
            typeof(h_pre),
            eltype(l_pre),
            typeof(res),
        }(em, h_pre, l_pre, res)
    end
end

function lscv(
    ::Type{D},
    X::AbstractVector,
    EM::EvaluationMethod=FFT();
    h0::Real = rule_of_thumb2(D,X),
    hlb::Real = h0*(1//10),
    hub::Real = h0*(5//2),
    hres::Integer = 20
) where {D<:Distribution}

    f(h::Real) = lscv(D, h, X, EM)
    h_pre = range(hlb, stop=hub, length=hres)
    l_pre = map(f, h_pre)

    ind = argmin(l_pre)
    if ind == 1
        error("""
            Hit bandwidth lower bound, try lowering hlb.
                Currently hlb=$hlb
            """
        )
    end
    if ind == hres
        error("""
            Hit bandwidth upper bound, try raising hub.
                Currently hub=$hub
            """
        )
    end

    h1 = h_pre[ind-1]
    h2 = h_pre[ind+1]

    res = optimize(f, h1, h2, GoldenSection(), store_trace=true)
    return UnivariateLSCV(D, EM, h_pre, l_pre, res)
end

@inline minimizer(ulscv::UnivariateLSCV) = Optim.minimizer(ulscv.res)
@inline Base.minimum(ulscv::UnivariateLSCV) = Optim.minimum(ulscv.res)

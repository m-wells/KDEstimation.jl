struct UnivariateLSCV{D,EM} <: LSCV{D,EM,1}
    eval_method::EM
    h_pre::StepRangeLen
    l_pre::Vector{Float64}
    res::UnivariateOptimizationResults{Float64,Float64,Float64,Float64,Float64,GoldenSection}

    UnivariateLSCV(::Type{D}, em::EM, h_pre::StepRangeLen, l_pre::Vector{Float64},
                   res::UnivariateOptimizationResults{Float64,Float64,Float64,Float64,Float64,GoldenSection}
                  ) where {D<:Distribution, EM<:EvaluationMethod} = new{D,EM}(em, h_pre, l_pre, res)
end

function lscv(::Type{D}, X::AbstractVector, EM::EvaluationMethod=FFT();
              h0::Float64 = rule_of_thumb2(D,X),
              hlb::Real = 0.1*h0,
              hub::Real = 2.5*h0,
              hres::Int = 20
             ) where {D<:Distribution}
    f(h::Float64) = lscv(D, h, X, EM)

    h_pre = range(hlb, stop=hub, length=hres)
    l_pre = f.(h_pre)

    ind = argmin(l_pre)
    ind == 1 && error("Hit lower boundary, try lowering hlb which is currently ", hlb)
    ind == hres && error("Hit upper boundary, try raising hub which is currently ", hub)

    h1 = h_pre[ind-1]
    h2 = h_pre[ind+1]

    res = optimize(f, h1, h2, GoldenSection(), store_trace=true)
    return UnivariateLSCV(D, EM, h_pre, l_pre, res)
end

@inline minimizer(ulscv::UnivariateLSCV) = Optim.minimizer(ulscv.res)
@inline Base.minimum(ulscv::UnivariateLSCV) = Optim.minimum(ulscv.res)

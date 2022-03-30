module KDEstimation

using StatsBase: quantile
using Distributions
using FFTW
using FFTW: plan_rfft, rFFTWPlan
using AbstractFFTs: ScaledPlan
using Interpolations
using Interpolations: FilledExtrapolation
using Optim

# the following line is needed for Tricube (or any other "custom" distribution)
using Distributions: @check_args, @distr_support
# can remove the above using line once distributions have been pulled into Distributions.jl

export kde, Direct, Binned, FFT
export rule_of_thumb1, rule_of_thumb2, lscv, minimizer, δ0, canonical_transform
export Biweight, Cosine, Epanechnikov, Normal, Triangular, Triweight, Uniform

include("univariate/kernels/kernels.jl")
include("univariate/kernels/Biweight.jl")
include("univariate/kernels/Cosine.jl")
include("univariate/kernels/Epanechnikov.jl")
include("./univariate/kernels/Logistic.jl")
include("univariate/kernels/Normal.jl")
include("univariate/kernels/Triangular.jl")
#include("./univariate/kernels/Tricube.jl")
include("univariate/kernels/Triweight.jl")
include("univariate/kernels/Uniform.jl")

abstract type EvaluationMethod end
struct Direct <: EvaluationMethod end
struct Binned <: EvaluationMethod
    M::Int
    Binned(M::Int) = new(M)
    Binned() = Binned(4096)
end
struct FFT <: EvaluationMethod
    M::Int

    FFT(M) = ispow2(M) ? new(M) : error("M needs to be a power of 2")
    FFT() = FFT(4096)
end

abstract type KDE{D<:Distribution, EM<:EvaluationMethod, N} end

function Base.getproperty(kde::KDE{D,EM,N}, s::Symbol) where {D,EM,N}
     s === :Distribution && return D
     s === :EvaluationMethod && return EM
     s === :N && return N
     return getfield(kde, s)
end

struct DirectUnivariateKDE{D} <: KDE{D,Direct,1}
    Kh::Kernel{Scaled,D}
    x::Vector{Float64}
end

struct BinnedUnivariateKDE{D} <: KDE{D,Binned,1}
    Kh::Kernel{Scaled,D}
    interp::FilledExtrapolation
end

struct FFTUnivariateKDE{D} <: KDE{D,FFT,1}
    Kh::Kernel{Scaled,D}
    interp::FilledExtrapolation
end

#Base.show(io::IO, ::KDE{D,Direct,N}) where {D,N} = print(io, "KDE{", D, ",Direct,",N,"}")
#
#function Base.show(io::IO, kde::KDE{D,EM,N}) where {D,EM,N}
#    nbins = length.(kde.interp.itp.ranges)
#    print(io, "KDE{", D, ",", EM)
#    length(nbins) == 1 ? print(io, "(", first(nbins), ")") : print(io, nbins)
#    print(io, ",", N, "}")
#end

function Base.show(io::IO, kde::KDE)
    print(io, "KDE{", join((kde.Distribution, kde.EvaluationMethod, kde.N), ", "), "}")
end

abstract type LSCV{
    D<:Distribution,
    EM<:EvaluationMethod,
    Res<:Optim.OptimizationResults,
    N,
} end

function Base.show(io::IO, lscv::LSCV{D,EM,Res,N}) where {D,EM,Res,N}
    println(io, "LSCV{",D,",",lscv.eval_method,",",N,"}")
    print(io, lscv.res)
end

include("./univariate/direct_kde.jl")
include("./univariate/binned_kde.jl")
include("./univariate/fft_kde.jl")

include("./univariate/bandwidth_selectors/canonical.jl")
include("./univariate/bandwidth_selectors/rot.jl")
include("./univariate/bandwidth_selectors/lscv.jl")

# Default to FFT method
function kde(D, h::Real, X::AbstractVector, EM=FFT())
    K = Kernel(D, h)
    return kde(K, X, EM)
end

function lscv(D, h::Real, X::AbstractVector, EM=FFT())
    K = Kernel(D, h)
    return lscv(K, X, EM)
end

include("./plot_recipes.jl")

end

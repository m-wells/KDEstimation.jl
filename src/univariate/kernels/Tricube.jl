# TODO: submit Tricube to Distributions.jl
struct Tricube{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    σ::T
    Tricube{T}(μ::T, σ::T) where {T<:Real} = new{T}(μ, σ)
end

function Tricube(μ::T, σ::T) where {T<:Real}
    @check_args(Tricube, σ ≥ zero(σ))
    return Tricube{T}(μ,σ)
end

#### Outer constructors
Tricube(μ::T, σ::T, ::NoArgCheck) where {T<:Real} = Tricube{T}(μ, σ)
Tricube(μ::Real, σ::Real) = Tricube(promote(μ,σ)...)

@distr_support Tricube (d.μ - d.σ) (d.μ + d.σ)

function pdf(d::Tricube, x::Real)
    u = abs(x - d.μ) / d.σ
    u ≥ 1 ? zero(T) : (70//81) * (1-u^3)^3 / d.σ
end#

# everything above here would go into Distributions.jl



@inline function autoconv(::Type{Tricube}, u::Real)
	x = abs(u)
	return x ≥ 2 ? 0.0 : 
           x ≥ 1 ? -((1/606092058)*(35*(-2 + x)^7*(151470     + 206822*x   + 558880*x^2 +
                                                   453722*x^3 + 561120*x^4 + 334740*x^5 +
                                                   213558*x^6 + 98448*x^7  + 33474*x^8  +
                                                   8439*x^9   + 1568*x^10  + 196*x^11   +
                                                   14*x^12))
                    ) :
                    (1/606092058)*(35*(12269070     - 19446804*x^2 + 23279256*x^4 - 51802740*x^6 +
                                       69006366*x^7 - 42854994*x^8 + 14965236*x^9 - 2863718*x^10 +
                                       71706*x^13   - 969*x^16     + 42*x^19)
                                  )
end

@inline μ2(::Type{Tricube}) = 35//243

@inline R(::Type{Tricube}) = 175//247

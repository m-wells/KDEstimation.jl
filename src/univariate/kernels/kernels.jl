abstract type KernelType end
struct Scaled <:KernelType end
struct UnScaled <:KernelType end

struct Kernel{S<:KernelType,D<:Distribution}
    d::D
    τ::Float64
    Kernel{S,D}(d::D,τ::Real) where {S<:KernelType, D<:Distribution} = new{S,D}(d,float(τ))
end

function Base.show(io::IO, k::K) where {K<:Kernel}
    println(io, K, "(τ => ", k.τ, ")")
end

"""
    get_τ(::Type{D}) where D<:Distribution = one(Float64)

τ provides the effective support.
For distributions that have infinite support (like the Normal) will need to supply
a more specific function.
"""
@inline get_τ(::Type{D}) where D<:Distribution = one(Float64)

@inline Kernel(::Type{D}) where D<:Distribution =
    Kernel{UnScaled,D}(D(0,1), get_τ(D))

@inline Kernel(::Type{D}, h::Real) where D<:Distribution =
    Kernel{Scaled,D}(D(0,h), get_τ(D)*h)

#@inline Scale(K::Kernel{UnScaled,D}, h::Real) where D<:Distribution = Kernel{Scaled,D}(D(0,h), K.τ*h)

@inline (K::Kernel)(u::Real) = pdf(K.d, u)

#@inline bandwidth(K::Kernel{S,D}) where {S<:KernelType, D<:Distribution} = K.d.σ
@inline bandwidth(K::Kernel{S,D}) where {S<:KernelType, D<:Distribution} = Distributions.scale(K.d)

autoconv(K::Kernel{UnScaled,D}, u::Real) where D<:Distribution = autoconv(D,u)

function autoconv(Kh::Kernel{Scaled,D}, u::Real) where D<:Distribution
    h = bandwidth(Kh)
    return autoconv(D,u/h)/h
end

function T_H(Kh::Kernel{Scaled,D}, u::Real) where {D<:Distribution}
    return autoconv(Kh,u) - 2*Kh(u)
end

@inline μ2(K::Kernel{S,D}) where {S<:KernelType, D<:Distribution} = μ2(D)
@inline R(K::Kernel{S,D}) where {S<:KernelType, D<:Distribution} = R(D)

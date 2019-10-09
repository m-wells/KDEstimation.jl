@inline function autoconv(::Type{Cosine}, u::Real)
    x = abs(u)
    return x ≥ 2 ? 0.0 : (-π/32)*(π*(-2 + x)*cos((π*x)/2) - 2*sin((π*x)/2))
end

@inline μ2(::Type{Cosine}) = 1 - 8/(π^2)

@inline R(::Type{Cosine}) = (π^2)/16

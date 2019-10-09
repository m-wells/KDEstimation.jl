@inline function autoconv(::Type{Epanechnikov}, u::Real)
    x = abs(u)
    return x ≥ 2 ? 0.0 : (-3/160)*(-2 + x)^3*(4 + 6*x + x^2)
end

@inline μ2(::Type{Epanechnikov}) = 1/5

@inline R(::Type{Epanechnikov}) = 3/5

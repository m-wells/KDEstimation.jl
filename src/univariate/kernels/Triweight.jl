@inline function autoconv(::Type{Triweight}, u::Real)
    x = abs(u)
    return x ≥ 2 ? 0.0 : -((35*(-2 + x)^7*(320 + 1120*x + 1616*x^2 + 1176*x^3 + 404*x^4 + 70*x^5 + 5*x^6))/1757184)
end

@inline μ2(::Type{Triweight}) = 1/9

@inline R(::Type{Triweight}) = 350/429

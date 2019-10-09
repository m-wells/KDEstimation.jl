const Triangular = SymTriangularDist

@inline function autoconv(::Type{Triangular}, u::Real)
    x = abs(u)
    return x ≥ 2 ? 0.0 : x ≥ 1 ? -(1/6)*(-2 + x)^3 : (1/6)*(4 - 6*x^2 + 3*x^3)
end

@inline μ2(::Type{Triangular}) = 1/6

@inline R(::Type{Triangular}) = 2/3

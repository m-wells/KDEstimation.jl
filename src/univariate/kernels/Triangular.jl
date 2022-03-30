const Triangular = SymTriangularDist

function autoconv(::Type{Triangular}, u::Real)
    x = abs(u)
    return if x ≥ 2
        0.0
    elseif x ≥ 1
        -(1/6)*(-2 + x)^3
    else
        (1/6)*(4 - 6*x^2 + 3*x^3)
    end
end

μ2(::Type{Triangular}) = 1/6

R(::Type{Triangular}) = 2/3

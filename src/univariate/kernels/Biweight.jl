#TODO submit scale(Biweight) upstream
@inline Distributions.scale(D::Biweight) = D.σ

@inline function autoconv(::Type{Biweight}, u::Real)
    x = abs(u)
    return x ≥ 2 ? 0.0 : -(5*(-2 + x)^5*(16 + 40*x + 36*x^2 + 10*x^3 + x^4))/3584
end

@inline μ2(::Type{Biweight}) = 1/7

@inline R(::Type{Biweight}) = 5/7

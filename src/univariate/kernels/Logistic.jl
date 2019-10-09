@inline get_τ(::Type{Logistic}, tol=1e-4) = abs(log((1 - sqrt(1-4*tol) - 2*tol)/(2*tol)))

@inline function autoconv(::Type{Logistic}, u::Real)
    # (exp(u)*(2 + exp(u)*(-2 + u) + u))/(-1 + exp(u))^3
    x = abs(u)

    # because of the exponentials roundoff becomes an issue
    # above -11.4572 the value of the autoconv drops below 10^-4
    x > 11.4572 && return 0.0

    # need to use L'Hôpital's rule because denominator vanishes at 0
    # using third derivatives
    numer = exp(x)*(5 + x + exp(x)*(-4 + 8x))
    denom = 3exp(x)*(1 - 8exp(x) + 9exp(2x))
    return numer/denom
end

@inline μ2(::Type{Logistic}) = (π^2)/3

@inline R(::Type{Logistic}) = 1/6

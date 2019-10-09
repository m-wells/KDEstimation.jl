@inline get_τ(::Type{Normal}, tol=1e-4) = sqrt(-2*log(tol*sqrt(2π)))

@inline autoconv(::Type{Normal}, u::Real) = pdf(Normal(0,√2), u)

@inline μ2(::Type{Normal}) = 1

@inline R(::Type{Normal}) = 1/(2*√π)

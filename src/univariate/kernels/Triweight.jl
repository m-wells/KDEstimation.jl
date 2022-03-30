function autoconv(::Type{Triweight}, u::T) where T
    x = abs(u)
    return if x ≥ 2
        zero(T)
    else
        #-((35*(-2 + x)^7*(320 + 1120*x + 1616*x^2 + 1176*x^3 + 404*x^4 + 70*x^5 + 5*x^6))/1757184)
        -(35*(-2+x)^7)* evalpoly(x, (320, 1120, 1616, 1176, 404, 70, 5))/1757184
    end
end

@inline μ2(::Type{Triweight}) = 1/9

@inline R(::Type{Triweight}) = 350/429

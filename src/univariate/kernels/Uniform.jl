@inline Kernel(::Type{Uniform}) = Kernel{UnScaled,Uniform}(Uniform(-1,1),1)
@inline Kernel(::Type{Uniform}, h::Real) = Kernel{Scaled,Uniform}(Uniform(-h,h),h)

#@inline Scale(::Kernel{UnScaled,Uniform}, h::Real) = Kernel{Scaled,Uniform}(Uniform(-h,h),h)

@inline bandwidth(K::Kernel{S,Uniform}) where {S<:KernelType} = (K.d.b - K.d.a)/2

@inline function autoconv(::Type{Uniform}, u::Real)
    x = abs(u)
    return x ≥ 2 ? 0.0 : (2-x)/4
end

@inline μ2(::Type{Uniform}) = 1/3

@inline R(::Type{Uniform}) = 1/2

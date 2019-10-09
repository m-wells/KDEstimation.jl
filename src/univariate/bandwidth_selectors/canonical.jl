"""
    canonical bandwidth scaling factor
"""
@inline δ0(::Type{D}) where {D<:Distribution} = (R(D)/μ2(D)^2)^(1/5)

"""
    canonical bandwidth transform
"""
function canonical_transform(h1::Real, ::Type{D1}, ::Type{D2}) where {D1<:Distribution, D2<:Distribution}
    δ0_1 = δ0(D1)
    δ0_2 = δ0(D2)
    return (h1/δ0_1)*δ0_2
end

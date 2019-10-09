@inline (kde::BinnedUnivariateKDE)(x) = kde.interp(x)

@inline function get_grid(Kh::Kernel, X::AbstractVector, M::Int)
    X1,X2 = extrema(X)
    return range(X1-Kh.τ, stop=X2+Kh.τ, length=M)
end

@inline function get_L(Kh::Kernel, grid::AbstractRange)
    M,δ = length(grid),step(grid)
    return min(M-1, ceil(Int, Kh.τ/δ))
end

"""
    get_2L(Kh::Kernel, grid::AbstractRange)

used for computations involving the autoconv
"""
@inline function get_2L(Kh::Kernel, grid::AbstractRange)
    M,δ = length(grid),step(grid)
    return min(M-1, 2*ceil(Int, Kh.τ/δ))
end

@inline function linear_binning(grid::AbstractRange, X::AbstractVector)
    M,δ = length(grid), step(grid)
    c = zeros(M)
    
    # linear binning
    @simd for x in X
        i = searchsortedlast(grid,x)
    
        if (i == 1) || (i == M)
            c[i] += 1
        else
            g1 = grid[i]
            g2 = grid[i+1]
    
            c[i] += (g2 - x)/δ
            c[i+1] += (x - g1)/δ
        end
    end
    return c
end

"""
Equation 5.19 from Nonparametric Kernel Density Estimation and Its Computational Aspects
"""
@inline function ckjl(c::AbstractVector, k::AbstractVector, j::Int, l::Int, L::Int)
    jl = j-l
    c_jl = jl in eachindex(c) ? c[jl] : 0.0
    k_l = k[L+1+l] # k_-L => k[1], k_0 => k[L+1], k_L => k[2L+1]
    return c_jl*k_l
end

function kde(Kh::Kernel{Scaled,D}, X::AbstractVector, EM::Binned) where {D<:Distribution}
    M = EM.M
    grid = get_grid(Kh,X,M)
    δ = step(grid)
    L = get_L(Kh, grid)

    c = linear_binning(grid,X)
    n = length(X)
    k = [Kh(δ*l)/n for l in -L:L]

    dgrid = Vector{Float64}(undef, M)
    for j in eachindex(grid)
        dgrid[j] = sum(l -> ckjl(c,k,j,l,L), -L:L)
    end

    return BinnedUnivariateKDE{D}(Kh, LinearInterpolation(grid,dgrid,extrapolation_bc=0))
end

function lscv(Kh::Kernel{Scaled,D}, X::AbstractVector, EM::Binned) where {D<:Distribution}
    M = EM.M
    grid = get_grid(Kh,X,M)
    δ = step(grid)
    L = get_2L(Kh, grid)

    c = linear_binning(grid,X)
    n = length(X)
    k = [T_H(Kh,δ*l) for l in -L:L]

    dgrid = Vector{Float64}(undef, M)
    for j in eachindex(grid)
        dgrid[j] = sum(l -> ckjl(c,k,j,l,L), -L:L)
    end

    ψ = sum(c.*dgrid)/(n^2)

    return ψ + 2*Kh(0)/n
end

lscv(::Type{D}, h::Real, X::AbstractVector, EM::Binned) where {D<:Distribution} = lscv(Kernel(D,h), X, EM)

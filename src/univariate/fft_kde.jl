const FFTPlansDict = Dict{Int64,Tuple{rFFTWPlan{Float64,-1,false,1},
                                      ScaledPlan{Complex{Float64},rFFTWPlan{Complex{Float64},1,false,1},Float64}
                                     }
                         }

global fftplans = FFTPlansDict()

function get_plans!(plans::FFTPlansDict, N::Int)
    N in keys(plans) && return plans[N]

    FFT = plan_rfft(Vector{Float64}(undef,N), flags=FFTW.MEASURE)
    IFFT = inv(FFT)
    plans[N] = (FFT,IFFT)
    return (FFT,IFFT)
end

@inline (kde::FFTUnivariateKDE)(x) = kde.interp(x)

@inline function kde(Kh::Kernel{Scaled,D}, X::AbstractVector, EM::FFT) where {D<:Distribution}
    M = EM.M
    grid = get_grid(Kh,X,M)
    δ = step(grid)
    L = get_L(Kh, grid)
    P = 2^ceil(Int, log2(M+2L+1))

    FFT,IFFT = get_plans!(fftplans::FFTPlansDict, P)

    cpz = zeros(P)
    cpz[(L+1):(L+M)] .= linear_binning(grid, X)

    n = length(X)

    kpz = zeros(P)
    for (j,l) in enumerate(-L:L)
        kpz[j] = Kh(δ*l)/n
    end

    C = FFT*cpz
    K = FFT*kpz
    C.*=K

    s = IFFT*C
    d_estimate = s[(2L+1):(2L+M)]

    return FFTUnivariateKDE(Kh, LinearInterpolation(grid, d_estimate, extrapolation_bc = 0))
end


@inline function lscv(Kh::Kernel{Scaled,D}, X::AbstractVector, EM::FFT) where {D<:Distribution}
    M = EM.M
    grid = get_grid(Kh,X,M)
    δ = step(grid)
    L = get_2L(Kh, grid)
    P = 2^ceil(Int, log2(M+2L+1))

    FFT,IFFT = get_plans!(fftplans::FFTPlansDict, P)

    c = linear_binning(grid, X)

    cpz = zeros(P)
    cpz[(L+1):(L+M)] .= c

    n = length(X)

    kpz = zeros(P)
    for (j,l) in enumerate(-L:L)
        kpz[j] = T_H(Kh,δ*l)
    end

    C = FFT*cpz
    C .*= FFT*kpz
    #K = FFT*kpz
    #C.*=K

    s = IFFT*C
    dgrid = s[(2L+1):(2L+M)]
    ψ = sum(c.*dgrid)/(n^2)

    return ψ + 2Kh(0)/n
end

@inline lscv(::Type{D}, h::Real, X::AbstractVector, EM::FFT) where {D<:Distribution} = lscv(Kernel(D,h), X, EM)

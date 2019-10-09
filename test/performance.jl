using KDEstimation
using Distributions
using Random: seed!


function direct_kde(x,t)
    dkde = kde(Normal,0.2,x,Direct())
    return dkde(t)
end

function binned_kde(x,t)
    bkde = kde(Normal,0.2,x,Binned())
    return bkde(t)
end

function fft_kde(x,t)
    fkde = kde(Normal,0.2,x,FFT())
    return fkde(t)
end


function direct_lscv(x,hs)
    return [lscv(Normal,h,x,Direct()) for h in hs]
end

function binned_lscv(x,hs)
    return [lscv(Normal,h,x,Binned()) for h in hs]
end

function fft_lscv(x,hs)
    return [lscv(Normal,h,x,FFT()) for h in hs]
end


function kde_times_compare()
    seed!(1234)
    x = randn(1000)
    t = range(-3,stop=3,length=1000)

    println("direct: ", direct_kde(x,0.1))
    println("binned: ", binned_kde(x,0.1))
    println("fft:    ", fft_kde(x,0.1))

    println("sample size = ", length(x))
    println("evaluations = ", length(t))
    println("direct")
    @time direct_kde(x,t)
    println("binned")
    @time binned_kde(x,t)
    println("fft")
    @time fft_kde(x,t)
    return nothing
end

function lscv_times_compare()
    seed!(1234)
    x = randn(1000)
    h = range(0.1,stop=3,length=100)

    @show direct_lscv(x,0.2)
    @show binned_lscv(x,0.2)
    @show fft_lscv(x,0.2)

    println("direct")
    @time direct_lscv(x,h)
    println("binned")
    @time binned_lscv(x,h)
    println("fft")
    @time fft_lscv(x,h)
    return nothing
end

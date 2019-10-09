using RecipesBase

@recipe function f(kde::DirectUnivariateKDE)
    h = bandwidth(kde.Kh)
    @series begin
        x1,x2 = extrema(kde.X)
        pad = 1.1*h
        x1 -= pad
        x2 += pad
        x = range(x1, stop=x2, length=1000)
        y = kde(x)
        seriestype --> :line
        x, y
    end
    @series begin
        seriestype --> :scatter
        markershape --> :vline
        kde.X, zeros(length(kde.X))
    end
end

@recipe function f(kde::K) where K<:Union{BinnedUnivariateKDE,FFTUnivariateKDE}
    h = bandwidth(kde.Kh)
    X = first(kde.interp.itp.ranges)
    Δx = step(X)
    @series begin
        x1,x2 = extrema(X)
        pad = 1.1*h
        x = [range(x1-pad, stop=x1-Δx, length=10); X; range(x2+Δx, stop=x2+pad, length=10)]
        y = kde(x)
        seriestype := :line
        x, y
    end
    #@series begin
    #    seriestype := :sticks
    #    X, Y
    #end
end



@recipe function f(ulscv::UnivariateLSCV{D,M}) where {D<:Distribution, M<:EvaluationMethod}
    n = ulscv.res.f_calls

    htrace = Vector{Float64}(undef,n)
    ltrace = Vector{Float64}(undef,n)
    for (i,t) in enumerate(ulscv.res.trace)
        htrace[i] = t.metadata["minimizer"]
        ltrace[i] = t.value
    end

    @series begin
        seriestype := :line
        label := string("lscv (", M, ")")
        [htrace; ulscv.h_pre], [ltrace; ulscv.l_pre]
    end

    @series begin
        seriestype := :scatter
        label := "presearch"
        ulscv.h_pre, ulscv.l_pre
    end

    @series begin
        seriestype := :scatter
        label := "solver trace"
        htrace,ltrace
    end

    hlb,hub = extrema(ulscv.h_pre)
    @series begin
        seriestype := :vline
        label := "bounds"
        [hlb,hub]
    end

    hmin = Optim.minimizer(ulscv.res)
    @series begin
        seriestype := :vline
        label := "optimial"
        [hmin]
    end

    lmin = Optim.minimum(ulscv.res)
    @series begin
        seriestype := :hline
        label := "minimum"
        [lmin]
    end
end

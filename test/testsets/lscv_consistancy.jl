dists = [# method, opt bwidth, opt lscv, tolerance
    (Biweight    , 1.5360    ,  -0.1320, 1e-3),
    (Cosine      , 3.9140    ,  -0.1385, 1e-3),
    (Epanechnikov, 1.1310    ,  -0.1322, 1e-3),
    (Logistic    , 0.3123    ,  -0.1683, 1e-3),
    (Normal      , 0.6987    ,  -0.1318, 1e-2),
    (Triangular  , 1.4580    ,  -0.1318, 1e-2),
    (Triweight   , 1.7480    ,  -0.1318, 1e-3),
    (Uniform     , 0.7690    ,  -0.1356, 1e-2),
]

ms = [Direct(), Binned(), FFT()]

x = 2 .* randn(StableRNG(123), 100) .- 5

@testset "$d" for (d,h,l,t) in dists
    for m in ms
        mlscv = lscv(d,x,m)
        #println(
        #    rpad(d, 20),
        #    "h = ",
        #    rpad(Float16(minimizer(mlscv)), 8),
        #    "l = ",
        #    Float16(minimum(mlscv)),
        #)
        #break
        @test isapprox(minimizer(mlscv), h, atol=t)
        @test isapprox(minimum(mlscv), l, atol=t)
    end
end

dists = [(Biweight    , 0.2109, 1e-3),
         (Cosine      , 0.2127, 1e-3),
         (Epanechnikov, 0.2041, 1e-3),
         (Logistic    , 0.1998, 1e-3),
         (Normal      , 0.2013, 1e-3),
         (Triangular  , 0.2083, 1e-3),
         (Triweight   , 0.2155, 1e-3),
         (Uniform     , 0.1969, 1e-2)
        ]

ms = [Direct(), Binned(), FFT()]

x = 2 .* randn(StableRNG(123), 1000) .- 5
h = 0.1

@testset "$d" for (d,v,tol) in dists
    for m in ms
        mkde = kde(d, h, x, m)
        #println(
        #    rpad(d, 20),
        #    "kde = ",
        #    mkde(-4.5),
        #)
        #break
        @test isapprox(mkde(-4.5), v, atol=tol)
    end
end

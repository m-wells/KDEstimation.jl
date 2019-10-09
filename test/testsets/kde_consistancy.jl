dists = [(Biweight    , 0.200, 1e-3),
         (Cosine      , 0.199, 1e-3),
         (Epanechnikov, 0.201, 1e-3),
         (Logistic    , 0.200, 1e-3),
         (Normal      , 0.203, 1e-3),
         (Triangular  , 0.203, 1e-3),
         (Triweight   , 0.199, 1e-3),
         (Uniform     , 0.197, 3e-3)
        ]

ms = [Direct(), Binned(), FFT()]

seed!(1234)
x = 2 .* randn(1000) .- 5
h = 0.1

@testset "$d" for (d,v,tol) in dists
    for m in ms
        mkde = kde(d, h, x, m)
        @test isapprox(mkde(-4.5), v, atol=tol)
    end
end

dists = [(Biweight    , 2.0362),
         (Cosine      , 1.7663),
         (Epanechnikov, 1.7188),
         (Logistic    , 0.4340),
         (Normal      , 0.7764),
         (Triangular  , 1.8882), # (Tricube     , 2.0262),
         (Triweight   , 2.3122),
         (Uniform     , 1.3510)
        ]

@testset "$d" for (d,v) in dists
    @test isapprox(δ0(d), v, atol=1e-4)
end

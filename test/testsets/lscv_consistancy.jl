seed!(1234)
x = 2 .* randn(100) .- 5

# method, optim bandwidth, lscv at optim, tolerance
dists = [(Biweight     , 2.556 , -0.1458 , 1e-3),
         (Cosine       , 3.261 , -0.1535 , 1e-3),
         (Epanechnikov , 2.015 , -0.1461 , 1e-3),
         (Logistic     , 0.257 , -0.1844 , 1e-3),
         (Normal       , 0.971 , -0.1453 , 1e-3),
         (Triangular   , 2.310 , -0.1459 , 1e-3),
         (Triweight    , 2.904 , -0.1457 , 1e-3),
         (Uniform      , 1.496 , -0.1495 , 1e-3)]

ms = [Direct(), Binned(), FFT()]

@testset "$d" for (d,h,l,t) in dists
    for m in ms
        mlscv = lscv(d,x,m)
        #@show d
        #@show minimizer(mlscv)
        #@show minimum(mlscv)
        @test isapprox(minimizer(mlscv), h, atol=t)
        @test isapprox(minimum(mlscv), l, atol=t)
    end
end

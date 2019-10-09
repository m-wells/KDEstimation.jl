using KDEstimation
using Test

using Random: seed!
using Distributions


@testset "canonical bandwidth" begin
    include("./testsets/cannon_bandwidths.jl")
end

@testset "kde consistancy" begin
    include("./testsets/kde_consistancy.jl")
end

@testset "lscv consistancy" begin
    include("./testsets/lscv_consistancy.jl")
end



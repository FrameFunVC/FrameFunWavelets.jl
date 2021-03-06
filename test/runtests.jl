using Test

@testset "DyadicPeriodicEquispacedGrids" begin
    include("test_dyadicgrids.jl")
end

@testset "WaveletBases" begin
    include("test_waveletbases.jl")
end

@testset "WaveletBasesPlatforms" begin
    include("test_waveletbasesplatforms.jl")
end

@testset "WaveletPlatforms" begin
    include("test_waveletplatforms.jl")
end

@testset "Compact AZ" begin
    include("test_compactaz.jl")
end

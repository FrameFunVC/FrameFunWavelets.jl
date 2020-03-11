module TestWaveletPlatforms

using Test, FrameFunWavelets, FrameFun, DomainSets, LinearAlgebra
@testset "plunge rank and nonzero columns" begin
    nnzcolumn(P,N,threshold) =
        sum(sum(abs.(plungematrix(P,N));dims=1).>threshold)

    nnzrow(P,N,threshold) =
        sum(sum(abs.(plungematrix(P,N));dims=2).>threshold)
    myplungerank(P,N,threshold) =
        sum(svdvals(plungematrix(P,N)) .> threshold)

    P = ExtensionFramePlatform(DaubechiesPlatform(2),0.0..0.4)
    @test myplungerank(P,3,1e-10) == 3
    @test myplungerank(P,4,1e-10) == 3
    @test myplungerank(P,5,1e-10) == 2
    @test myplungerank(P,6,1e-10) == 3
    @test myplungerank(P,7,1e-10) == 3

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 13
    @test nnzcolumn(P,5,1e-10) == 17
    @test nnzcolumn(P,6,1e-10) == 23
    @test nnzcolumn(P,7,1e-10) == 26

    @test nnzrow(P,3,1e-10) == 17
    @test nnzrow(P,4,1e-10) == 19
    @test nnzrow(P,5,1e-10) == 7
    @test nnzrow(P,6,1e-10) == 20
    @test nnzrow(P,7,1e-10) == 17

    P = ExtensionFramePlatform(DaubechiesPlatform(3),0.0..0.4)
    @test myplungerank(P,3,1e-10) == 5
    @test myplungerank(P,4,1e-10) == 5
    @test myplungerank(P,5,1e-10) == 5
    @test myplungerank(P,6,1e-10) == 5
    @test myplungerank(P,7,1e-10) == 5

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 16
    @test nnzcolumn(P,5,1e-10) == 25
    @test nnzcolumn(P,6,1e-10) == 33
    @test nnzcolumn(P,7,1e-10) == 42

    @test nnzrow(P,3,1e-10) == 26
    @test nnzrow(P,4,1e-10) == 35
    @test nnzrow(P,5,1e-10) == 19
    @test nnzrow(P,6,1e-10) == 36
    @test nnzrow(P,7,1e-10) == 33

    # Is it numeric noise that increases plunge rank?
    P = ExtensionFramePlatform(DaubechiesPlatform(5),0.0..0.4)
    @test myplungerank(P,3,1e-10) == 8
    @test myplungerank(P,4,1e-10) == 13
    @test myplungerank(P,5,1e-10) == 14
    @test myplungerank(P,6,1e-10) == 13
    @test myplungerank(P,7,1e-10) == 13

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 16
    @test nnzcolumn(P,5,1e-10) == 31
    @test nnzcolumn(P,6,1e-10) == 48
    @test nnzcolumn(P,7,1e-10) == 65

    @test nnzrow(P,3,1e-10) == 26
    @test nnzrow(P,4,1e-10) == 52
    @test nnzrow(P,5,1e-10) == 52
    @test nnzrow(P,6,1e-10) == 99
    @test nnzrow(P,7,1e-10) == 96

    P = ExtensionFramePlatform(CDFPlatform(2,2),0.0..0.4)
    @test myplungerank(P,3,1e-10) == 1
    @test myplungerank(P,4,1e-10) == 1
    @test myplungerank(P,5,1e-10) == 1
    @test myplungerank(P,6,1e-10) == 1
    @test myplungerank(P,7,1e-10) == 1

    @test nnzcolumn(P,3,1e-10) == 6
    @test nnzcolumn(P,4,1e-10) == 9
    @test nnzcolumn(P,5,1e-10) == 13
    @test nnzcolumn(P,6,1e-10) == 15
    @test nnzcolumn(P,7,1e-10) == 17

    @test nnzrow(P,3,1e-10) == 1
    @test nnzrow(P,4,1e-10) == 2
    @test nnzrow(P,5,1e-10) == 4
    @test nnzrow(P,6,1e-10) == 3
    @test nnzrow(P,7,1e-10) == 1

    P = ExtensionFramePlatform(CDFPlatform(3,3),0.0..0.4)
    @test myplungerank(P,3,1e-10) == 3
    @test myplungerank(P,4,1e-10) == 3
    @test myplungerank(P,5,1e-10) == 4
    @test myplungerank(P,6,1e-10) == 4
    @test myplungerank(P,7,1e-10) == 3

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 16
    @test nnzcolumn(P,5,1e-10) == 26
    @test nnzcolumn(P,6,1e-10) == 36
    @test nnzcolumn(P,7,1e-10) == 44

    @test nnzrow(P,3,1e-10) == 12
    @test nnzrow(P,4,1e-10) == 13
    @test nnzrow(P,5,1e-10) == 15
    @test nnzrow(P,6,1e-10) == 14
    @test nnzrow(P,7,1e-10) == 12

    P = ExtensionFramePlatform(CDFPlatform(6,6),0.0..0.4)
    @test myplungerank(P,3,1e-8) == 6
    @test myplungerank(P,4,1e-8) == 7
    @test myplungerank(P,5,1e-8) == 6
    @test myplungerank(P,6,1e-8) == 7
    @test myplungerank(P,7,1e-8) == 7

    @test nnzcolumn(P,3,1e-8) == 8
    @test nnzcolumn(P,4,1e-8) == 16
    @test nnzcolumn(P,5,1e-8) == 32
    @test nnzcolumn(P,6,1e-8) == 53
    @test nnzcolumn(P,7,1e-8) == 74

    @test nnzrow(P,3,1e-8) == 17
    @test nnzrow(P,4,1e-8) == 32
    @test nnzrow(P,5,1e-8) == 29
    @test nnzrow(P,6,1e-8) == 33
    @test nnzrow(P,7,1e-8) == 31
end



end

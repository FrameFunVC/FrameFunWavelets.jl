module TestWaveletPlatforms

using Test, FrameFunWavelets, FrameFun, DomainSets
@testset "AZ_Zt" begin
    P = DaubechiesPlatform(2)
    N = 4
    op = AZ_Z(P,N)
    @test element(op,1) isa InverseDiscreteWaveletTransform
    @test element(op,2) isa VerticalBandedOperator
    op2 = AZ_Zt(P,N)

    @test element(op2,length(elements(op))) isa DiscreteWaveletTransform
    @test element(op2,length(elements(op))-1) isa HorizontalBandedOperator
    @test AZ_Zt(P,N)*AZ_A(P,N) ≈ IdentityOperator(dictionary(P,N))

    P = CDFPlatform(4,4)
    N = 4
    g = sampling_grid(P,N)
    dual_dict = dualdictionary(P,N,discretemeasure(g))
    dict = dictionary(P,N)
    g = sampling_grid(P,N)
    Z = AZ_Z(P,N)
    E = evaluation_operator(dual_dict, g)
    @test Z≈E
    A = AZ_A(P,N)
    E = evaluation_operator(dict, g)
    @test A≈E

    @test Z'A ≈ IdentityOperator(dict)
    @test AZ_Zt(P,N)*AZ_A(P,N) ≈ IdentityOperator(dict)
end


@testset "plunge rank and nonzero columns" begin

    nnzcolumn(P,N,threshold) =
        sum(sum(abs.(plungematrix(P,N));dims=1).>threshold)

    P = ExtensionFramePlatform(DaubechiesPlatform(2),0.0..0.4)
    @test plungerank(P,3;threshold=1e-10) == 3
    @test plungerank(P,4;threshold=1e-10) == 3
    @test plungerank(P,5;threshold=1e-10) == 2
    @test plungerank(P,6;threshold=1e-10) == 3
    @test plungerank(P,7;threshold=1e-10) == 3

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 13
    @test nnzcolumn(P,5,1e-10) == 17
    @test nnzcolumn(P,6,1e-10) == 23
    @test nnzcolumn(P,7,1e-10) == 26

    P = ExtensionFramePlatform(DaubechiesPlatform(3),0.0..0.4)
    @test plungerank(P,3;threshold=1e-10) == 5
    @test plungerank(P,4;threshold=1e-10) == 5
    @test plungerank(P,5;threshold=1e-10) == 5
    @test plungerank(P,6;threshold=1e-10) == 5
    @test plungerank(P,7;threshold=1e-10) == 5

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 16
    @test nnzcolumn(P,5,1e-10) == 25
    @test nnzcolumn(P,6,1e-10) == 33
    @test nnzcolumn(P,7,1e-10) == 40

    P = ExtensionFramePlatform(DaubechiesPlatform(5),0.0..0.4)
    @test plungerank(P,3;threshold=1e-8) == 5
    @test plungerank(P,4;threshold=1e-8) == 9
    @test plungerank(P,5;threshold=1e-8) == 10
    @test plungerank(P,6;threshold=1e-8) == 9
    @test plungerank(P,7;threshold=1e-8) == 9

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 16
    @test nnzcolumn(P,5,1e-10) == 31
    @test nnzcolumn(P,6,1e-10) == 52
    @test nnzcolumn(P,7,1e-10) == 86

    P = ExtensionFramePlatform(CDFPlatform(2,2),0.0..0.4)
    @test plungerank(P,3;threshold=1e-10) == 1
    @test plungerank(P,4;threshold=1e-10) == 1
    @test plungerank(P,5;threshold=1e-10) == 1
    @test plungerank(P,6;threshold=1e-10) == 1
    @test plungerank(P,7;threshold=1e-10) == 1

    @test nnzcolumn(P,3,1e-10) == 7
    @test nnzcolumn(P,4,1e-10) == 9
    @test nnzcolumn(P,5,1e-10) == 11
    @test nnzcolumn(P,6,1e-10) == 15
    @test nnzcolumn(P,7,1e-10) == 19

    P = ExtensionFramePlatform(CDFPlatform(3,3),0.0..0.4)
    @test plungerank(P,3;threshold=1e-10) == 3
    @test plungerank(P,4;threshold=1e-10) == 3
    @test plungerank(P,5;threshold=1e-10) == 2
    @test plungerank(P,6;threshold=1e-10) == 3
    @test plungerank(P,7;threshold=1e-10) == 3

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 16
    @test nnzcolumn(P,5,1e-10) == 24
    @test nnzcolumn(P,6,1e-10) == 35
    @test nnzcolumn(P,7,1e-10) == 45

    P = ExtensionFramePlatform(CDFPlatform(6,6),0.0..0.4)
    @test plungerank(P,3;threshold=1e-8) == 6
    @test plungerank(P,4;threshold=1e-8) == 7
    @test plungerank(P,5;threshold=1e-8) == 6
    @test plungerank(P,6;threshold=1e-8) == 6
    @test plungerank(P,7;threshold=1e-8) == 7

    @test nnzcolumn(P,3,1e-10) == 8
    @test nnzcolumn(P,4,1e-10) == 16
    @test nnzcolumn(P,5,1e-10) == 32
    @test nnzcolumn(P,6,1e-10) == 54
    @test nnzcolumn(P,7,1e-10) == 76

end


end

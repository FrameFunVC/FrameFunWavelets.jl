module TestWaveletPlatforms
using Test, FrameFun, FrameFunWavelets, DomainSets, LinearAlgebra
using WaveletsEvaluation.DWT: CDFWavelet, DaubechiesWavelet, Dul, Wvl, Prl

@testset "Constructors" begin
    @test platform(DaubechiesWaveletBasis(3,10)) isa DaubechiesPlatform{3,Float64,Wvl}
    @test platform(CDFWaveletBasis(3,2,10)) isa CDFPlatform{3,2,Float64,Prl,Wvl}
    @test CDFPlatform{Float64}(3,2) isa CDFPlatform{3,2,Float64,Prl,Wvl}
    @test DaubechiesPlatform{Float64}(3) isa DaubechiesPlatform{3,Float64,Wvl}
    @test CDFPlatform(3,2) isa CDFPlatform{3,2,Float64,Prl,Wvl}
    @test DaubechiesPlatform(3) isa DaubechiesPlatform{3,Float64,Wvl}
    @test Platform(CDFWavelet{3,2,Float64}()) isa CDFPlatform{3,2,Float64,Prl,Wvl}
    @test Platform(DaubechiesWavelet{3,Float64}()) isa DaubechiesPlatform{3,Float64,Wvl}
end

@testset "dictionary" begin
    P = DaubechiesPlatform(3)
    @test dictionary(P,4) isa DaubechiesWaveletBasis
    @test azdual_dict(P,4;samplingstyle=GramStyle()) isa DaubechiesWaveletBasis

    P = CDFPlatform(4,2)
    @test dictionary(P,4) isa CDFWaveletBasis
    @test azdual_dict(P,4;samplingstyle=GramStyle()) isa CDFWaveletBasis
end

@testset "AZ_A" begin
    P = DaubechiesPlatform(3)
    @test dictionary(P,10) isa DaubechiesWaveletBasis
    @test sampling_grid(P,10) == interpolation_grid(dictionary(P,10))
    @test dualdictionary(P,10,FourierMeasure()) == dictionary(P,10)
    @test element(AZ_A(P,5),1) isa DWTEvalOperator

    PP = ExtensionFramePlatform(P,0.0..0.49)
    g = sampling_grid(PP,5)
    @test supergrid(g) isa DyadicPeriodicEquispacedGrid
    @test element(AZ_A(PP,5), 1) isa DWTEvalOperator

    P = CDFPlatform(4,2)
    @test dictionary(P,10) isa CDFWaveletBasis{4,2,Float64,Prl,Wvl}
    @test sampling_grid(P,10) == interpolation_grid(dictionary(P,10))
    @test dualdictionary(P,10,FourierMeasure()) isa CDFWaveletBasis{4,2,Float64,Dul,Wvl}
    @test element(AZ_A(P,5),1) isa InverseDiscreteWaveletTransform

    PP = ExtensionFramePlatform(P,0.0..0.49)
    g = sampling_grid(PP,5)
    @test supergrid(g) isa PeriodicEquispacedGrid
    @test length(supergrid(g))/32==4
    g = sampling_grid(PP,5;oversamplingfactor=3)
    @test supergrid(g) isa PeriodicEquispacedGrid
    @test length(supergrid(g))/32==6
    @test element(AZ_A(P,5),1) isa InverseDiscreteWaveletTransform
    op = AZ_A(PP,5)
    @test element(op, 1) isa InverseDiscreteWaveletTransform
    @test element(op, 2) isa VerticalBandedOperator
end

@testset "Direct approximation" begin
    f = exp; N = 5
    P = ExtensionFramePlatform(CDFPlatform(2,2),(0.0..0.5)^1)
    F = Fun(f, P, N; solverstyle=DirectStyle())
    @test norm(AZ_A(P,N)*coefficients(F)-f.(sampling_grid(P,N))) < 1e-3

    f = (x,y)->exp(x*y); N = (5,5)
    P = ExtensionFramePlatform(NdCDFPlatform(2,2,2),(0.0..0.5)^2)
    F = Fun(f, P, N; solverstyle=DirectStyle())
    @test norm(AZ_A(P,N)*coefficients(F)-[f(x...) for x in sampling_grid(P,N)]) < 1e-3

    f = (x,y,z)->exp(x*y*z); N = (3,3,3)
    P = ExtensionFramePlatform(NdCDFPlatform(3,2,2),(0.0..0.5)^3)
    F = Fun(f, P, N; solverstyle=DirectStyle())
    @test norm(AZ_A(P,N)*coefficients(F)-[f(x...) for x in sampling_grid(P,N)]) < 1e-3

    f = exp; N = 5
    P = ExtensionFramePlatform(DaubechiesPlatform(2),(0.0..0.5)^1)
    F = Fun(f, P, N; solverstyle=DirectStyle())
    @test norm(AZ_A(P,N)*coefficients(F)-f.(sampling_grid(P,N))) < 1e-2


    f = (x,y)->exp(x*y); N = (5,5)
    P = ExtensionFramePlatform(NdDaubechiesPlatform(2,2),(0.0..0.5)^2)
    F = Fun(f, P, N; solverstyle=DirectStyle())
    @test norm(AZ_A(P,N)*coefficients(F)-[f(x...) for x in sampling_grid(P,N)]) < 1e-2

    f = (x,y,z)->exp(x*y*z); N = (3,3,3)
    P = ExtensionFramePlatform(NdDaubechiesPlatform(3,2),(0.0..0.5)^3)
    F = Fun(f, P, N; solverstyle=DirectStyle())
    @test norm(AZ_A(P,N)*coefficients(F)-[f(x...) for x in sampling_grid(P,N)]) < 1e-2
end


using FrameFunWavelets.WaveletBasesPlatforms.CompactPeriodicEquispacedWaveletDual: signal
@testset "signal" begin
    N = 4
    P = CDFPlatform(2,2)
    dict = dictionary(P,N)
    @test evaluation_operator(scalingbasis(dict), PeriodicEquispacedGrid(2*length(dict),support(dict))).A[:,1]==signal(scalingbasis(dict), 2)[0:length(dict)*2-1]

    N = 4
    P = DaubechiesPlatform(2)
    dict = dictionary(P,N)
    g = PeriodicEquispacedGrid(2*length(dict),support(dict))
    @test evaluation_operator(scalingbasis(dict), g).A[:,1]≈signal(scalingbasis(dict), 2)[0:length(dict)*2-1]
end

using FrameFunWavelets, FrameFun, Test, LinearAlgebra, InfiniteVectors
using BSplineExtension.BSplinePlatforms.CompactPeriodicEquispacedTranslatesDuals: CompactPeriodicEquispacedTranslatesDual, signal
using FrameFunWavelets.WaveletPlatforms.WaveletBasesPlatforms.CompactPeriodicEquispacedWaveletDual: compact_wavelet_dual
@testset "dual dictionary" begin
    N = 4
    P = DaubechiesPlatform(2)
    dict = scalingbasis(dictionary(P,N))
    g = DyadicPeriodicEquispacedGrid(2(1<<N),FrameFun.support(dict))
    E1 = evaluation_operator(dict, g).A
    b = signal(dict, 2)
    c = PeriodicInfiniteVector(b,32)[0:31]
    d = PeriodicInfiniteVector(inv(b, 2;K=1)',32)[0:31]
    @test dot(c,d)≈1

    N = 4
    P = DaubechiesPlatform(2)
    dict = dictionary(P,N)
    g = DyadicPeriodicEquispacedGrid(2(1<<N),FrameFun.support(dict))
    dual_dict = dualdictionary(P,N,discretemeasure(g))

    E1 = evaluation_operator(dict, g)
    S1 = evaluation_operator(scalingbasis(dict), g)
    iDWT1 = InverseDiscreteWaveletTransform(dict)
    @test S1*iDWT1≈E1

    E2 = evaluation_operator(dual_dict, PeriodicEquispacedGrid(g))
    S2 = element(E2,2)
    iDWT2 = element(E2,1)

    @test Matrix(S1.A'S2.A)≈I
    @test Matrix(iDWT1)'*Matrix(iDWT2)≈I
    @test Matrix(E1)'Matrix(E2)≈I

    N = 4
    P = CDFPlatform(4,4)
    dict = dictionary(P,N)
    g = DyadicPeriodicEquispacedGrid(2(1<<N),FrameFun.support(dict))
    dual_dict = dualdictionary(P,N,discretemeasure(g))

    E1 = evaluation_operator(dict, g)
    S1 = evaluation_operator(scalingbasis(dict), g)
    iDWT1 = InverseDiscreteWaveletTransform(dict)
    @test S1*iDWT1≈E1

    E2 = evaluation_operator(dual_dict, PeriodicEquispacedGrid(g))
    S2 = element(E2,2)
    iDWT2 = element(E2,1)

    @test Matrix(S1.A'S2.A)≈I
    @test Matrix(iDWT1)'*Matrix(iDWT2)≈I
    @test Matrix(E1)'Matrix(E2)≈I
end


end
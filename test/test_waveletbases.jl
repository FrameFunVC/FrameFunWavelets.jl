module TestWaveletBases
using Test, FrameFun.BasisFunctions, FrameFun.BasisFunctions.Test, FrameFunWavelets.WaveletBases,
    StaticArrays, DomainSets, FrameFunWavelets.DyadicPeriodicEquispacedGrids

using WaveletsEvaluation.DWT: quad_trap, quad_sf, quad_sf_weights, quad_sf_N, quad_trap_N
using WaveletsEvaluation.DWT: wavelet, Dual, scaling, db3,  db4, Primal, Prl, value

@testset "DWT, iDWT" begin
    for dict in (CDFWaveletBasis(4,2,5),CDFWaveletBasis(4,4,10),DaubechiesWaveletBasis(3,5),
                    CDFWaveletBasis(4,2,5,Prl,Float64,false),CDFWaveletBasis(4,4,10,Prl,Float64,false),DaubechiesWaveletBasis(3,5,Float64,false))

        DWT1 = DiscreteWaveletTransform(dict)
        iDWT1 = InverseDiscreteWaveletTransform(dict)
        DWT2 = DiscreteWaveletTransform(wavelet_dual(dict))
        iDWT2 = InverseDiscreteWaveletTransform(wavelet_dual(dict))


        @test inv(iDWT1)≈DWT1
        @test inv(DWT1)≈iDWT1
        @test inv(iDWT2)≈DWT2
        @test inv(DWT2)≈iDWT2
        @test iDWT1'≈DWT2
        @test iDWT2'≈DWT1
        @test DWT1'≈iDWT2
        @test DWT2'≈iDWT1

        @test inv(Matrix(iDWT1))≈Matrix(DWT1)
        @test inv(Matrix(DWT1))≈Matrix(iDWT1)
        @test inv(Matrix(iDWT2))≈Matrix(DWT2)
        @test inv(Matrix(DWT2))≈Matrix(iDWT2)
        @test Matrix(iDWT1)'≈Matrix(DWT2)
        @test Matrix(iDWT2)'≈Matrix(DWT1)
        @test Matrix(DWT1)'≈Matrix(iDWT2)
        @test Matrix(DWT2)'≈Matrix(iDWT1)
    end
end

@testset "wavelet quadrature" begin
    M1 = 5; M2 = 10; J = 3
    wav = db3; g = x->sin(2pi*x)
    reference = quad_sf_N(g, wav, M2, J, 5)

    b = DaubechiesWaveletBasis(3,J)
    S = DWTSamplingOperator(b, Int(M1/5),0)
    @test norm(S*g-reference/sqrt(1<<J)) < 1e-3
    S = DWTSamplingOperator(b, Int(M2/5), 0)
    @test norm(S*g-reference/sqrt(1<<J)) < 1e-8
    S = DWTSamplingOperator(b, Int(M1/5),7)
    @test norm(S*g-reference/sqrt(1<<J)) < 1e-15
    S = DWTSamplingOperator(b, Int(M2/5), 3)
    @test norm(S*g-reference/sqrt(1<<J)) < 1e-15
end

@testset "DaubechiesWaveletBasis and CDFWaveletBasis" begin

    b1 = DaubechiesWaveletBasis(3,2)
    b2 = CDFWaveletBasis(3,1,5)
    b = CDFWaveletBasis(1,1,3)
    BasisFunctions.unsafe_eval_element(b, 1, .1)
    160 == @allocated BasisFunctions.unsafe_eval_element(b, 1, .1)
    supports = ((0.,1.),(0.,1.),(0.0,0.5),(0.5,1.0),(0.0,0.25),(0.25,0.5),(0.5,0.75),(0.75,1.0));
    for i in ordering(b)
        @test infimum(support(b,i)) == supports[value(i)][1]
        @test supremum(support(b,i)) == supports[value(i)][2]
    end
    for i in ordering(b1)
        @test infimum(support(b1,i))==0
        @test supremum(support(b1,i))==1
    end
    @test BasisFunctions.subdict(b1,1:1) == DaubechiesWaveletBasis(3,0)
    @test BasisFunctions.subdict(b2,1:3) == BasisFunctions.LargeSubdict(b2,1:3)
    @test b2[1:5] == BasisFunctions.LargeSubdict(b2,1:5)
    @test b2[1:4] == CDFWaveletBasis(3,1,2)
    @test dyadic_length(b1) == 2
    @test dyadic_length(b2) == 5
    @test length(b1) == 4
    @test length(b2) == 32
    @test BasisFunctions.promote_domaintype(b1,ComplexF64) == DaubechiesWaveletBasis(3,2, ComplexF64)
    @test BasisFunctions.promote_domaintype(b2, ComplexF64) == CDFWaveletBasis(3,1,5, Prl, ComplexF64)
    @test resize(b1,8) == DaubechiesWaveletBasis(3,3)
    @test BasisFunctions.name(b1) == "Basis of db3 wavelets"
    @test BasisFunctions.name(b2) == "Basis of cdf31 wavelets"
    @test BasisFunctions.native_index(b,1) == (scaling, 0, 0)
    @test BasisFunctions.native_index(b,2) == (wavelet, 0, 0)
    @test BasisFunctions.native_index(b,3) == (wavelet, 1, 0)
    @test BasisFunctions.native_index(b,4) == (wavelet, 1, 1)
    @test BasisFunctions.native_index(b,5) == (wavelet, 2, 0)
    @test BasisFunctions.native_index(b,6) == (wavelet, 2, 1)
    @test BasisFunctions.native_index(b,7) == (wavelet, 2, 2)
    @test BasisFunctions.native_index(b,8) == (wavelet, 2, 3)
    for i in 1:length(b)
        @test(linear_index(b,BasisFunctions.native_index(b,i))==i )
    end
    @test interpolation_grid(b1) == PeriodicEquispacedGrid(4,0,1)
    @test interpolation_grid(b2) == PeriodicEquispacedGrid(32,0,1)
    @test BasisFunctions.period(b1)==1.

    # test grid eval functions
    for g in (plotgrid(b,128), PeriodicEquispacedGrid(128,0,1))
        for i in ordering(b)
            t=@timed e1 = BasisFunctions._default_unsafe_eval_element_in_grid(b, i, g)
            t1 = t[2]
            t = @timed e2 = WaveletBases._unsafe_eval_element_in_dyadic_grid(b, i, g)
            t2 = t[2]
            @test e1 ≈ e2
        end
        for i in ordering(b)
            t = @timed e1 = BasisFunctions._default_unsafe_eval_element_in_grid(b, i, g)
            t1 = t[2]
            t = @timed e2 = WaveletBases._unsafe_eval_element_in_dyadic_grid(b, i, g)
            t2 = t[2]
            @test t2 < t1
        end
    end

    b1 = DaubechiesWaveletBasis(3,6)
    b2 = CDFWaveletBasis(3,1,5)
    b = CDFWaveletBasis(1,1,6)
    for basis in (b,b1,b2, wavelet_dual(b),scalingbasis(b))
        @test evaluation_matrix(basis, interpolation_grid(basis))≈evaluation_matrix(basis, collect(interpolation_grid(basis)))
        # T = Float64
        # ELT = Float64
        # tbasis = transform_dict(basis)
        # x = ones(length(basis))
        # t = transform_operator(tbasis, basis)
        # it = transform_operator(basis, tbasis)
        # @test norm((it*t*x-x))+1≈1
end

dict = DaubechiesWaveletBasis(2,10)
g = interpolation_grid(dict)
@test g isa DyadicPeriodicEquispacedGrid
op = evaluation_operator(dict, g)
@test element(op,1) isa InverseDiscreteWaveletTransform
op = evaluation_operator(dict, PeriodicEquispacedGrid(interpolation_grid(dict)))
@test element(op,1) isa InverseDiscreteWaveletTransform

dict = CDFWaveletBasis(2,2,10)
g = interpolation_grid(dict)
@test g isa DyadicPeriodicEquispacedGrid
op = evaluation_operator(dict, g)
@test element(op,1) isa InverseDiscreteWaveletTransform
@test element(op,2) isa VerticalBandedOperator
op = evaluation_operator(dict, PeriodicEquispacedGrid(interpolation_grid(dict)))
@test element(op,1) isa InverseDiscreteWaveletTransform
@test element(op,2) isa VerticalBandedOperator
end

@testset "Scaled wavelet basis" begin
dict = DaubechiesWaveletBasis(2,5,Float64,false)
DWT = DiscreteWaveletTransform(dict)
iDWT = InverseDiscreteWaveletTransform(dict)

inv(DWT)≈DWT
Matrix(iDWT*iDWT')
inv(iDWT)


@test DWT*iDWT≈IdentityOperator(dict)
dict2 = DaubechiesWaveletBasis(2,5,Float64,true)
DWT2 = DiscreteWaveletTransform(dict2)
iDWT2 = InverseDiscreteWaveletTransform(dict2)
@test DWT2*iDWT2≈IdentityOperator(dict2)
@test !(DWT≈DWT2)
S1 = evaluation_operator(scalingbasis(dict),interpolation_grid(dict))
S2 = evaluation_operator(scalingbasis(dict2),interpolation_grid(dict))
@test S1.A.array*sqrt(32)≈S2.A.array
W1 = evaluation_operator(dict,interpolation_grid(dict))
W2 = evaluation_operator(dict2,interpolation_grid(dict))
e = zeros(dict)
e[1] = 1
@test W1*e≈W2*e
@test W1 ≈ S1*iDWT
@test W2 ≈ S2*iDWT2

dict = CDFWaveletBasis(2,4,5,Prl,Float64,false)
DWT = DiscreteWaveletTransform(dict)
iDWT = InverseDiscreteWaveletTransform(dict)
@test DWT*iDWT≈IdentityOperator(dict)
dict2 = CDFWaveletBasis(2,4,5,Prl,Float64,true)
DWT2 = DiscreteWaveletTransform(dict2)
iDWT2 = InverseDiscreteWaveletTransform(dict2)
@test DWT2*iDWT2≈IdentityOperator(dict2)
@test !(DWT≈DWT2)
S1 = evaluation_operator(scalingbasis(dict),interpolation_grid(dict))
S2 = evaluation_operator(scalingbasis(dict2),interpolation_grid(dict))
@test S1.A.array*sqrt(32)≈S2.A.array
W1 = evaluation_operator(dict,interpolation_grid(dict))
W2 = evaluation_operator(dict2,interpolation_grid(dict))
e = zeros(dict)
e[1] = 1
@test W1*e≈W2*e
@test W1 ≈ S1*iDWT
@test W2 ≈ S2*iDWT2
end

@testset "Plots" begin
    using Plots
    b = CDFWaveletBasis(1,5,3)
    plot(b,layout=2)
    plot(wavelet_dual(b)[3:4],layout=2,subplot=1)
    plot!(Dual, wavelet, wavelet(b),j=1,k=0,subplot=2,periodic=true)
    plot!(Dual, wavelet, wavelet(b),j=1,k=1,subplot=2,periodic=true)
end
end

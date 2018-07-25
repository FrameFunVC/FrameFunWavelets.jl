# test_wavelets.jl
using BasisFunctions, Domains, StaticArrays, WaveletsDict, BasisFunctions.Test
if VERSION < v"0.7-"
    using Base.Test
    ComplexF64 = Complex128
else
    using Test
end

using WaveletsCopy.DWT: quad_trap, quad_sf, quad_sf_weights, quad_sf_N, quad_trap_N
using WaveletsCopy.DWT: wavelet, Dual, scaling, db3,  db4, Primal, Prl, value

for T in (Float64, BigFloat)
    @testset "$(rpad(string(PeriodicEquispacedGrid),80))" begin
        g = instantiate(PeriodicEquispacedGrid,10,T)
        test_interval_grid(g)
    end
end


function test_wavelet_quadrature()
        @testset begin
            M1 = 5; M2 = 10; J = 3
            wav = db3; g = x->sin(2pi*x)
            reference = quad_sf_N(g, wav, M2, J, 5)

            b = DaubechiesWaveletBasis(3,J)
            S = WaveletsDict.DWTSamplingOperator(b, Int(M1/5),0)
            @test norm(S*g-reference/sqrt(1<<J)) < 1e-3
            S = WaveletsDict.DWTSamplingOperator(b, Int(M2/5), 0)
            @test norm(S*g-reference/sqrt(1<<J)) < 1e-8
            S = WaveletsDict.DWTSamplingOperator(b, Int(M1/5),7)
            @test norm(S*g-reference/sqrt(1<<J)) < 1e-15
            S = WaveletsDict.DWTSamplingOperator(b, Int(M2/5), 3)
            @test norm(S*g-reference/sqrt(1<<J)) < 1e-15
        end
end

test_wavelet_quadrature()


function bf_wavelets_implementation_test()
    @testset begin
        # test_generic_dict_interface(CDFWaveletBasis(3,5,6))
        # test_generic_dict_interface(wavelet_dual(CDFWaveletBasis(3,5,6)))
        # test_generic_dict_interface(CDFScalingBasis(1,1,6))

        b1 = DaubechiesWaveletBasis(3,2)
        b2 = CDFWaveletBasis(3,1,5)
        b = CDFWaveletBasis(1,1,3)
        BasisFunctions.unsafe_eval_element(b, 1, .1)
        624 == @allocated BasisFunctions.unsafe_eval_element(b, 1, .1)
        supports = ((0.,1.),(0.,1.),(0.0,0.5),(0.5,1.0),(0.0,0.25),(0.25,0.5),(0.5,0.75),(0.75,1.0));
        for i in ordering(b)
            @test infimum(support(b,i)) == supports[value(i)][1]
            @test supremum(support(b,i)) == supports[value(i)][2]
        end
        for i in ordering(b1)
            @test support(b1,i) == UnitInterval()
        end
        @test BasisFunctions.subdict(b1,1:1) == DaubechiesWaveletBasis(3,0)
        @test BasisFunctions.subdict(b2,1:3) == BasisFunctions.LargeSubdict(b2,1:3)
        @test b2[1:5] == BasisFunctions.LargeSubdict(b2,1:5)
        @test b2[1:4] == CDFWaveletBasis(3,1,2)
        @test dyadic_length(b1) == 2
        @test dyadic_length(b2) == 5
        @test length(b1) == 4
        @test length(b2) == 32
        @test BasisFunctions.dict_promote_domaintype(b1,ComplexF64) == DaubechiesWaveletBasis(3,2, ComplexF64)
        @test BasisFunctions.dict_promote_domaintype(b2, ComplexF64) == CDFWaveletBasis(3,1,5, Prl, ComplexF64)
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
        @test grid(b1) == PeriodicEquispacedGrid(4,0,1)
        @test grid(b2) == PeriodicEquispacedGrid(32,0,1)
        @test BasisFunctions.period(b1)==1.

        # test grid eval functions
        for g in (plotgrid(b,200), PeriodicEquispacedGrid(128,0,1))
            for i in ordering(b)
                t=@timed e1 = BasisFunctions._default_unsafe_eval_element_in_grid(b, i, g)
                t1 = t[2]
                t = @timed e2 = WaveletsDict._unsafe_eval_element_in_dyadic_grid(b, i, g)
                t2 = t[2]
            end
            for i in ordering(b)
                t = @timed e1 = BasisFunctions._default_unsafe_eval_element_in_grid(b, i, g)
                t1 = t[2]
                t = @timed e2 = WaveletsDict._unsafe_eval_element_in_dyadic_grid(b, i, g)
                t2 = t[2]
                @test e1 ≈ e2
                @test t2 < t1
            end
        end

        b1 = DaubechiesWaveletBasis(3,6)
        b2 = CDFWaveletBasis(3,1,5)
        b = CDFWaveletBasis(1,1,6)
        for basis in (b,b1,b2,wavelet_dual(b), ScalingBasis(b))
            @test evaluation_matrix(basis, BasisFunctions.grid(basis))≈evaluation_matrix(basis, collect(BasisFunctions.grid(basis)))
            T = Float64
            ELT = Float64
            tbasis = transform_dict(basis)
            x = ones(length(basis))
            t = transform_operator(tbasis, basis)
            it = transform_operator(basis, tbasis)
            @test norm((it*t*x-x))+1≈1
        end
    end
end
bf_wavelets_implementation_test()

if VERSION < v"0.7-"
    using Plots
    b = CDFWaveletBasis(1,5,3)
    plot(b,layout=2)
    plot(wavelet_dual(b)[3:4],layout=2,subplot=1)
    plot!(Dual, wavelet, wavelet(b),j=1,k=0,subplot=2,periodic=true)
    plot!(Dual, wavelet, wavelet(b),j=1,k=1,subplot=2,periodic=true)
end

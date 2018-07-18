using BasisFunctions, Domains, StaticArrays, WaveletsDict
using Base.Test
using WaveletsCopy: cdf33, db3, cdf13

# function test_coefficient_index_range_of_overlapping_elements()
#     L = 5
#     S = ScalingBasis(cdf33, L)
#     for x in [0, 0.01, 0.99]
#         g = ScatteredGrid([x])
#         @test sort(find(evaluation_matrix(S, g).!=0)) == sort(collect(coefficient_index_range_of_overlapping_elements(S, g[1])))
#     end
#     S = S⊗S
#     for x in [0, 0.01, 0.99]
#         g = ScatteredGrid([x])
#         g = g×g
#         @test sort(find(evaluation_matrix(S, g).!=0)) == sort([sub2ind(size(S), i.I...) for i in coefficient_index_range_of_overlapping_elements(S, g[1])])
#     end
# end


function test_scaling_platform()
    platform = scaling_platform([4,5], [db3,db3], 2)
    B = primal(platform, 1)
    e = rand(B)
    @test B == ScalingBasis(db3, 4)⊗ScalingBasis(db3, 5)
    @test dual(platform, 1) == wavelet_dual(B)
    @test sampler(platform, 1)==GridSamplingOperator( BasisFunctions.oversampled_grid(B, 2))
    @test dual_sampler(platform, 1).sampler==sampler(platform, 1)
    e = rand(src(dual_sampler(platform, 1).weight))
    @test dual_sampler(platform, 1).weight*e≈WaveletsDict.WeightOperator(primal(platform, 1), [2,2], [0,0])*e
    @test BasisFunctions.Zt(platform, 1)*(sampler(platform, 1)*((x,y)->1.))≈ones(B)/sqrt(length(B))
    e = rand(src(BasisFunctions.A(platform, 1)))
    @test BasisFunctions.A(platform, 1)*e≈evaluation_operator(B, BasisFunctions.oversampled_grid(B, 2))*e

    platform = scaling_platform([4,5], [db3,cdf13], 2)
    B = primal(platform, 2)
    @test B == ScalingBasis(db3, 5)⊗ScalingBasis(cdf13, 6)
    @test dual(platform, 2) == wavelet_dual(B)
    @test sampler(platform, 2)==GridSamplingOperator( BasisFunctions.oversampled_grid(B, 2))
    @test dual_sampler(platform, 2).sampler==sampler(platform, 2)
    e = rand(src(dual_sampler(platform, 2).weight))
    @test dual_sampler(platform, 2).weight*e≈WaveletsDict.WeightOperator(primal(platform, 2), [2,2], [0,0])*e
    @test BasisFunctions.Zt(platform, 2)*(sampler(platform, 2)*((x,y)->1.))≈ones(B)/sqrt(length(B))
    e = rand(src(BasisFunctions.A(platform, 2)))
    @test BasisFunctions.A(platform, 2)*e≈evaluation_operator(B, BasisFunctions.oversampled_grid(B, 2))*e

    platform = scaling_platform([4,5], [db3,cdf13], 4)
    B = primal(platform, 2)
    @test B == ScalingBasis(db3, 5)⊗ScalingBasis(cdf13, 6)
    @test dual(platform, 2) == wavelet_dual(B)
    @test sampler(platform, 2)==GridSamplingOperator( BasisFunctions.oversampled_grid(B, 4))
    @test dual_sampler(platform, 2).sampler==sampler(platform, 2)
    e = rand(src(dual_sampler(platform, 2).weight))
    @test dual_sampler(platform, 2).weight*e≈WaveletsDict.WeightOperator(primal(platform, 2), [2,2], [1,1])*e
    @test BasisFunctions.Zt(platform, 2)*(sampler(platform, 2)*((x,y)->1.))≈ones(B)/sqrt(length(B))
    e = rand(src(BasisFunctions.A(platform, 2)))
    @test BasisFunctions.A(platform, 2)*e≈evaluation_operator(B, BasisFunctions.oversampled_grid(B, 4))*e

end
# @testset "Scaling util (1)" begin test_coefficient_index_range_of_overlapping_elements() end
@testset "Scaling platform" begin test_scaling_platform() end

module TestDyadicGrids

using Test, GridArrays, DomainSets,
    FrameFunWavelets.DyadicPeriodicEquispacedGrids,
    BasisFunctions.Test

using GridArrays: rescale, hasextension, extend
using GridArrays.Test: test_generic_grid

for T in (Float64, BigFloat)
    @testset "$(rpad(string(DyadicPeriodicEquispacedGrid),80))" begin
        g = DyadicPeriodicEquispacedGrid{T}(4, 0, 1)
        test_generic_grid(g)
        T = eltype(g)
        g1 = rescale(rescale(g, -T(10), T(3)), infimum(covering(g)), supremum(covering(g)))
        @test covering(g1) ≈ covering(g)
        g2 = resize(g, length(g)<<1)
        @test length(g2) == length(g)<<1


        g3 = resize(g1, length(g)<<1)
        @test length(g3) == length(g)<<1

        g4 = rescale(rescale(g2, -T(10), T(3)), infimum(covering(g2)), supremum(covering(g2)))
        @test covering(g4) ≈ covering(g2)

        if hasextension(g)
            @test extend(g,4) isa typeof(g)
        end
    end
end

end # module

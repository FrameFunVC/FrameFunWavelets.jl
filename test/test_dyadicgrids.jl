module TestDyadicGrids
using Test, GridArrays, FrameFunWavelets.DyadicPeriodicEquispacedGrids, FrameFun.BasisFunctions.Test, DomainSets
using GridArrays: rescale, hasextension, extend

for T in (Float64, BigFloat)
    @testset "$(rpad(string(DyadicPeriodicEquispacedGrid),80))" begin
        g = DyadicPeriodicEquispacedGrid{T}(4)
        test_generic_grid(g)
        T = eltype(g)
        g1 = rescale(rescale(g, -T(10), T(3)), infimum(support(g)), supremum(support(g)))
        @test support(g1) ≈ support(g)
        g2 = resize(g, length(g)<<1)
        @test length(g2) == length(g)<<1


        g3 = resize(g1, length(g)<<1)
        @test length(g3) == length(g)<<1

        g4 = rescale(rescale(g2, -T(10), T(3)), infimum(support(g2)), supremum(support(g2)))
        @test support(g4) ≈ support(g2)

        if hasextension(g)
            @test extend(g,4) isa typeof(g)
        end
    end
end
end

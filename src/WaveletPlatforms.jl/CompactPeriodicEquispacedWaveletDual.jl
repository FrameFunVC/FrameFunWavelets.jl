module CompactPeriodicEquispacedWaveletDual


using ....FrameFunWavelets.WaveletBases, FrameFun.BasisFunctions
using ....FrameFunWavelets.DyadicPeriodicEquispacedGrids

using FrameFun.BasisFunctions: unsafe_matrix
using ....FrameFunWavelets.WaveletBases: GenericPeriodicEquispacedTranslates, isdyadic, dyadic_length
using WaveletsEvaluation.DWT: Scl, Wvl, Prl, evaluate_in_dyadic_points, scaling
using InfiniteVectors: sublength, CompactInfiniteVector

import BSplineExtension.BSplinePlatforms.CompactPeriodicEquispacedTranslatesDuals: CompactPeriodicEquispacedTranslatesDual, signal
import FrameFun.BasisFunctions: grid_evaluation_operator

compact_wavelet_dual(dict::WaveletBasis{T,S,Scl}, m) where {T,S} =
    CompactPeriodicEquispacedTranslatesDual(dict, m)

function compact_wavelet_dual(dict::WaveletBasis{T,S,Wvl}, m) where {T,S}
    dual_scalingdict = compact_wavelet_dual(scalingbasis(dict), m)
    dual_waveletdict = wavelet_dual(dict)
    InverseDiscreteWaveletTransform(dual_waveletdict, dual_scalingdict, wavelet(dict), side(dual_waveletdict), dyadic_length(dict))*dual_waveletdict
end

function signal(dict::WaveletBasis{T,S,Scl}, m::Int) where {T,S}
    @assert isdyadic(m)
    f, x = evaluate_in_dyadic_points(side(dict),scaling, wavelet(dict), dyadic_length(dict), 0, Int(log2(m))+dyadic_length(dict);points=true)
    b = CompactInfiniteVector(f, findfirst(x.==0)-1)
end

function signal(dict::CDFWaveletBasis{P,Q,T,Prl,Scl}, m::Int) where {P,Q,T}
    A = evaluation_operator(dict, PeriodicEquispacedGrid(m*length(dict),support(dict)))
    @assert A isa VerticalBandedOperator
    a = convert(Vector{T}, unsafe_matrix(A).array)
    f = 1
    l = length(a)
    for i in 1:length(a)
        if !(a[i] + 1 ≈ 1)
            break
        end
        f += 1
    end
    for j in length(a):-1:1
        if !(a[j] + 1 ≈ 1)
            break
        end
        l -= 1
    end
    truncatedarray = a[f:l]
    os = mod(unsafe_matrix(A).offset+f-1, m*length(dict))
    if !(0 <= os + length(truncatedarray) <= m*length(dict))
        os -= length(grid)
    end
    CompactInfiniteVector(truncatedarray, os)
end

grid_evaluation_operator(dict::CompactPeriodicEquispacedTranslatesDual, gb::GridBasis, grid::DyadicPeriodicEquispacedGrid; options...) =
    grid_evaluation_operator(dict, gb, PeriodicEquispacedGrid(grid))

end

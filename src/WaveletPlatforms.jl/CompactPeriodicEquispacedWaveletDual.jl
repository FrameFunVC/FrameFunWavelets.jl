module CompactPeriodicEquispacedWaveletDual

using ....FrameFunWavelets.WaveletBases, BasisFunctions, CompactTranslatesDict,
    ....FrameFunWavelets.DyadicPeriodicEquispacedGrids

using WaveletsEvaluation.DWT: Scl, Wvl, Prl, evaluate_in_dyadic_points, scaling

import BasisFunctions: evaluation

compact_wavelet_dual(dict::WaveletBasis{T,S,Scl}, m; options...) where {T,S} =
    CompactPeriodicEquispacedTranslatesDual(dict, m; options...)

function compact_wavelet_dual(dict::WaveletBasis{T,S,Wvl,scaled}, m; options...) where {T,S,scaled}
    dual_scalingdict = compact_wavelet_dual(scalingbasis(dict), m; options...)
    dual_waveletdict = wavelet_dual(dict)
    InverseDiscreteWaveletTransform(dual_waveletdict, dual_scalingdict, wavelet(dict), side(dual_waveletdict), dyadic_length(dict), scaled)*dual_waveletdict
end

evaluation(::Type{T}, dict::CompactPeriodicEquispacedTranslatesDual, gb::GridBasis, grid::DyadicPeriodicEquispacedGrid; options...) where {T} =
    evaluation(T, dict, gb, PeriodicEquispacedGrid(grid); options...)

end

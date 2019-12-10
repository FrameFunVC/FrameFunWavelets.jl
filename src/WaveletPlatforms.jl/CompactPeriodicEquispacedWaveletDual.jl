module CompactPeriodicEquispacedWaveletDual

using ....FrameFunWavelets.WaveletBases, BasisFunctions, CompactTranslatesDict,
    ....FrameFunWavelets.DyadicPeriodicEquispacedGrids

using WaveletsEvaluation.DWT: Scl, Wvl, Prl, evaluate_in_dyadic_points, scaling

import BasisFunctions: grid_evaluation_operator

compact_wavelet_dual(dict::WaveletBasis{T,S,Scl}, m) where {T,S} =
    CompactPeriodicEquispacedTranslatesDual(dict, m)

function compact_wavelet_dual(dict::WaveletBasis{T,S,Wvl,scaled}, m) where {T,S,scaled}
    dual_scalingdict = compact_wavelet_dual(scalingbasis(dict), m)
    dual_waveletdict = wavelet_dual(dict)
    InverseDiscreteWaveletTransform(dual_waveletdict, dual_scalingdict, wavelet(dict), side(dual_waveletdict), dyadic_length(dict), scaled)*dual_waveletdict
end

grid_evaluation_operator(dict::CompactPeriodicEquispacedTranslatesDual, gb::GridBasis, grid::DyadicPeriodicEquispacedGrid; options...) =
    grid_evaluation_operator(dict, gb, PeriodicEquispacedGrid(grid))

end

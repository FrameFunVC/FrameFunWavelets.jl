module FrameFunWavelets

using Reexport

include("DyadicPeriodicEquispacedGrids.jl")
@reexport using .DyadicPeriodicEquispacedGrids

include("WaveletBases.jl")
@reexport using .WaveletBases

include("WaveletPlatforms.jl/WaveletPlatforms.jl")
@reexport using .WaveletPlatforms

end

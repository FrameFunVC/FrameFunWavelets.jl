module FrameFunWavelets

using Reexport

include("DyadicPeriodicEquispacedGrids.jl")
@reexport using .DyadicPeriodicEquispacedGrids

include("WaveletBases.jl")
@reexport using .WaveletBases

include("sparseiDWTEs.jl")
using .sparseiDWTEs

include("WaveletPlatforms.jl/WaveletPlatforms.jl")
@reexport using .WaveletPlatforms

include("CompactAZ.jl")
@reexport using .CompactAZ

end

module WaveletPlatforms

using Reexport

include("WaveletBasesPlatforms.jl")
@reexport using .WaveletBasesPlatforms

include("ndplatforms.jl")

end

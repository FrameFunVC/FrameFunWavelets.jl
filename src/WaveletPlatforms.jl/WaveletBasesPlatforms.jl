module WaveletBasesPlatforms

include("CompactPeriodicEquispacedWaveletDual.jl")
using .CompactPeriodicEquispacedWaveletDual: compact_wavelet_dual

using FrameFunWavelets, FrameFun.Platforms, FrameFun.BasisFunctions
using WaveletsEvaluation.DWT: Wvl, Prl, Kind, Side, CDFWavelet, DaubechiesWavelet, DiscreteWavelet
import FrameFun.Platforms: platform, dictionary, SolverStyle, measure, SamplingStyle,
    dualdictionary, correctparamformat, unsafe_dictionary, Platform
import FrameFun.FrameFunInterface: correct_sampling_parameter, SamplingStrategy, oversampling_grid

export AbstractWaveletPlatform
abstract type AbstractWaveletPlatform{T,K,scaled} <: BasisPlatform end

SolverStyle(p::AbstractWaveletPlatform, ::SamplingStyle) = DualStyle()
SamplingStyle(::AbstractWaveletPlatform{T,K,false}) where {T,K} = OversamplingStyle()

correct_sampling_parameter(::SamplingStrategy, ::AbstractWaveletPlatform, param, L; options...) = error()
correct_sampling_parameter(::SamplingStrategy, ::AbstractWaveletPlatform, param::Int, L::Int; options...) =
    (round(Int, 1<<round(Int,log2(L/param))) * param)
correctparamformat(::AbstractWaveletPlatform, ::Int) = true
dualdictionary(platform::AbstractWaveletPlatform, param, measure::FourierMeasure; options...) =
    wavelet_dual(unsafe_dictionary(platform, param))

export scalingplatform
scalingplatform(P::AbstractWaveletPlatform) = platform(scalingbasis(dictionary(P,1)))
scalingplatform(P::ProductPlatform) = ProductPlatform(map(scalingplatform, elements(P)))
export waveletplatform
waveletplatform(P::AbstractWaveletPlatform) = platform(waveletbasis(dictionary(P,1)))
waveletplatform(P::ProductPlatform) = ProductPlatform(map(waveletplatform, elements(P)))
waveletplatform(w::DiscreteWavelet, ::Type{S}=Prl, scaled::Bool=false) where {S} = platform(waveletbasis(w,0,S,scaled))

function dualdictionary(platform::AbstractWaveletPlatform, param, measure::UniformDiracComb;
        options...)
    dict = dictionary(platform, param)
    g = grid(measure)
    @assert support(dict) ≈ support(g)
    m = length(g) / length(dict)
    @assert round(Int,m) ≈ m
    m = round(Int, m)
    if m == 1
        @warn "No compact dual possible, try oversampling"
        return dual(dict, measure; options...)
    else
        if isperiodic(g)
            return compact_wavelet_dual(dict, m; options...)
        else
            error()
        end
    end
end

export CDFPlatform
struct CDFPlatform{P,Q,T,S<:Side,K<:Kind,scaled} <: AbstractWaveletPlatform{T,K,scaled}
end

CDFPlatform(P::Int, Q::Int, ::Type{S}=Prl, ::Type{K}=Wvl, scaled::Bool=false) where {S<:Side,K<:Kind} =
    CDFPlatform{P,Q,Float64,S,K,scaled}()
CDFPlatform{T}(P::Int, Q::Int, ::Type{S}=Prl, ::Type{K}=Wvl, scaled::Bool=false) where {T,S<:Side,K<:Kind} =
    CDFPlatform{P,Q,T,S,K,scaled}()
Platform(w::CDFWavelet{P,Q,T}, ::Type{S}=Prl, ::Type{K}=Wvl, scaled::Bool=false) where {P,Q,T,S<:Side,K<:Kind} =
    CDFPlatform{P,Q,T,S,K,scaled}()


unsafe_dictionary(platform::CDFPlatform{P,Q,T,S,K,scaled}, param::Int) where {P,Q,T,S<:Side,K<:Kind,scaled} =
    CDFWaveletBasis{P,Q,T,S,K,scaled}(CDFWavelet{P,Q,T}(),param)
platform(::CDFWaveletBasis{P,Q,T,S,K,scaled}) where {P,Q,T,S,K,scaled} =
    CDFPlatform{P,Q,T,S,K,scaled}()
oversampling_grid(dict::CDFWaveletBasis, L) = PeriodicEquispacedGrid(L, support(dict))
correct_sampling_parameter(::SamplingStrategy, ::CDFPlatform{P,Q,T,Prl,K}, param::Int, L::Int; options...)  where {P,Q,T,K} =
    (round(Int, L/(1<<param)) * (1<<param))

export DaubechiesPlatform
struct DaubechiesPlatform{P,T,K<:Kind,scaled} <: AbstractWaveletPlatform{T,K,scaled}
end

DaubechiesPlatform(P::Int, ::Type{K}=Wvl, scaled::Bool=true) where {K} =
    DaubechiesPlatform{P,Float64,K,scaled}()
DaubechiesPlatform{T}(P::Int, ::Type{K}=Wvl, scaled::Bool=true) where {T,K} =
    DaubechiesPlatform{P,T,K,scaled}()
Platform(w::DaubechiesWavelet{P,T},::Type{K}=Wvl, scaled::Bool=true) where {P,T,K} =
    DaubechiesPlatform{P,T,K,scaled}()

unsafe_dictionary(platform::DaubechiesPlatform{P,T,K,scaled}, param::Int) where {P,T,K,scaled} =
    DaubechiesWaveletBasis{P,T,Prl,K,scaled}(DaubechiesWavelet{P,T}(), param)
platform(::DaubechiesWaveletBasis{P,T,S,K,scaled}) where {P,T,S,K,scaled} =
    DaubechiesPlatform{P,T,K,scaled}()

end

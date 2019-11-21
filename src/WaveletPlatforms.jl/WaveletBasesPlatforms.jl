module WaveletBasesPlatforms

include("CompactPeriodicEquispacedWaveletDual.jl")
using .CompactPeriodicEquispacedWaveletDual: compact_wavelet_dual

using FrameFunWavelets, FrameFun.Platforms, FrameFun.BasisFunctions
using WaveletsEvaluation.DWT: Wvl, Prl, Kind, Side, CDFWavelet, DaubechiesWavelet
import FrameFun.Platforms: platform, dictionary, SolverStyle, measure, SamplingStyle,
    dualdictionary, correctparamformat, unsafe_dictionary, Platform
import FrameFun.FrameFunInterface: correct_sampling_parameter, SamplingStrategy, oversampling_grid


abstract type AbstractWaveletPlatform{T,S} <: BasisPlatform end

SolverStyle(p::AbstractWaveletPlatform, ::SamplingStyle) = DualStyle()
correct_sampling_parameter(::SamplingStrategy, ::AbstractWaveletPlatform, param, L; options...) = error()
correct_sampling_parameter(::SamplingStrategy, ::AbstractWaveletPlatform, param::Int, L::Int; options...) =
    (round(Int, 1<<round(Int,log2(L/param))) * param)
correctparamformat(::AbstractWaveletPlatform, ::Int) = true
dualdictionary(platform::AbstractWaveletPlatform, param, measure::FourierMeasure; options...) =
    wavelet_dual(unsafe_dictionary(platform, param))

function dualdictionary(platform::AbstractWaveletPlatform, param, measure::UniformDiracCombMeasure;
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
            return compact_wavelet_dual(dict, m)
        else
            error()
        end
    end
end

export CDFPlatform
struct CDFPlatform{P,Q,T,S<:Side,K<:Kind} <: AbstractWaveletPlatform{T,T}
end

CDFPlatform(P::Int, Q::Int, ::Type{S}=Prl, ::Type{K}=Wvl) where {S<:Side,K<:Kind} =
    CDFPlatform{P,Q,Float64,S,K}()
CDFPlatform{T}(P::Int, Q::Int, ::Type{S}=Prl, ::Type{K}=Wvl) where {T,S<:Side,K<:Kind} =
    CDFPlatform{P,Q,T,S,K}()
Platform(w::CDFWavelet{P,Q,T}, ::Type{S}=Prl, ::Type{K}=Wvl) where {P,Q,T,S<:Side,K<:Kind} =
    CDFPlatform{P,Q,T,S,K}()

unsafe_dictionary(platform::CDFPlatform{P,Q,T,S,K}, param::Int) where {P,Q,T,S<:Side,K<:Kind} =
    CDFWaveletBasis{P,Q,T,S,K}(CDFWavelet{P,Q,T}(),param)
platform(::CDFWaveletBasis{P,Q,T,S,K}) where {P,Q,T,S,K} =
    CDFPlatform{P,Q,T,S,K}()
oversampling_grid(dict::CDFWaveletBasis, L) = PeriodicEquispacedGrid(L, support(dict))
correct_sampling_parameter(::SamplingStrategy, ::CDFPlatform{P,Q,T,Prl,K}, param::Int, L::Int; options...)  where {P,Q,T,K} =
    (round(Int, L/(1<<param)) * (1<<param))

export DaubechiesPlatform
struct DaubechiesPlatform{P,T,K<:Kind} <: AbstractWaveletPlatform{T,T}
end

DaubechiesPlatform(P::Int, ::Type{K}=Wvl) where {K} =
    DaubechiesPlatform{P,Float64,K}()
DaubechiesPlatform{T}(P::Int, ::Type{K}=Wvl) where {T,K} =
    DaubechiesPlatform{P,T,K}()
Platform(w::DaubechiesWavelet{P,T},::Type{K}=Wvl) where {P,T,K} =
    DaubechiesPlatform{P,T,K}()

unsafe_dictionary(platform::DaubechiesPlatform{P,T,K}, param::Int) where {P,T,K} =
    DaubechiesWaveletBasis{P,T,Prl,K}(DaubechiesWavelet{P,T}(), param)
platform(::DaubechiesWaveletBasis{P,T,S,K}) where {P,T,S,K} =
    DaubechiesPlatform{P,T,K}()

end

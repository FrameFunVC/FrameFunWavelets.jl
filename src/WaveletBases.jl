module WaveletBases

using Reexport

include("DWTOperators.jl")
@reexport using .DWTOperators
import .DWTOperators: DiscreteWaveletTransform, InverseDiscreteWaveletTransform

using BasisFunctions, GridArrays, DomainSets, ..DyadicPeriodicEquispacedGrids,
    WaveletsEvaluation, RecipesBase, InfiniteVectors, CompactTranslatesDict

using WaveletsEvaluation.DWT: Side, Kind, Prl, Scl, Wvl, DiscreteWavelet, WaveletIndex,
    wavelet_index, scaling_indices, scaling_index,
    evaluate_periodic_in_dyadic_points!, dwt!, offset, dyadicpointsofrecursion,
    evaluate_periodic_scaling_basis_in_dyadic_points
using BasisFunctions: SamplingOperator
using CardinalBSplines: evaluate_BSpline

import Base: checkbounds, convert, promote_rule
import BasisFunctions: length, size, native_index, name, hasextension,
    approx_length, resize, hasinterpolationgrid, hastransform, iscompatible, hasgrid_transform,
    support, measure, period, interpolation_grid, ordering, linear_index, extensionsize,
    unsafe_eval_element, unsafe_eval_element1,
    apply!, evaluation, apply, dest, src_space, plotgrid,
    promote_domaintype, isbasis, isorthogonal, isorthonormal, isbiorthogonal, evaluation_matrix!,
    ArrayOperator, symbol, _default_unsafe_eval_element_in_grid,
    isinplace, dest, diag, apply_inplace!, isdiag, adjoint
import ..DyadicPeriodicEquispacedGrids: dyadic_length
import WaveletsEvaluation.DWT: value, wavelet, kind
import CompactTranslatesDict: GenericPeriodicEquispacedTranslates


export WaveletBasis
abstract type WaveletBasis{T,S,K,scaled} <: Dictionary1d{T,T} where {S <: Side,K<:Kind}
end

export WaveletTensorProductDict
const WaveletTensorProductDict = TensorProductDict{N,<:Tuple{Vararg{<:WaveletBasis}}} where N

native_index(dict::WaveletBasis, idx::WaveletIndex) = idx

checkbounds(::Type{Bool}, dict::WaveletBasis, i::WaveletIndex) =
    checkbounds(Bool, dict, linear_index(dict, i))

"""
    dyadic_length(dict::WaveletBasis)
The number of levels in the wavelet basis
"""
dyadic_length(dict::WaveletBasis) = dict.L

length(dict::WaveletBasis) = 1<<dyadic_length(dict)
size(dict::WaveletBasis) = (length(dict),)

export wavelet
"""
    wavelet(dict::WaveletBasis)

The wavelet type of the basis.
"""
WaveletsEvaluation.DWT.wavelet(dict::WaveletBasis) = dict.w

name(b::WaveletBasis) = "Basis of "*WaveletsEvaluation.DWT.name(wavelet(b))*" wavelets"
export side
side(::WaveletBasis{T,S}) where {T,S} = S()
export kind
kind(::WaveletBasis{T,S,K}) where {T,S,K} = K()

hasextension(b::WaveletBasis{T,S}) where{T,S} = true

approx_length(b::WaveletBasis, n::Int) = 1<<ceil(Int, log2(n))

resize(b::B, n::Int) where {B<:WaveletBasis} = B(wavelet(b),round(Int, log2(n)))

hasinterpolationgrid(::WaveletBasis) = true

hastransform(::WaveletBasis) = true

iscompatible(set::WaveletBasis, grid::PeriodicEquispacedGrid) =
    length(set)==length(grid)
iscompatible(set::WaveletBasis, grid::DyadicPeriodicEquispacedGrid) =
	length(set)==length(grid)
hasgrid_transform(b::WaveletBasis, gb, grid) = iscompatible(b, grid)

support(dict::WaveletBasis{T}) where T = UnitInterval{T}()
measure(dict::WaveletBasis{T}) where T = FourierMeasure{T}()
support(b::WaveletBasis, idx) = support(b, native_index(b, idx))

function support(dict::WaveletBasis{T}, idxn::WaveletIndex) where {T}
    l,r = WaveletsEvaluation.support(side(dict), kind(idxn), wavelet(dict), level(idxn), offset(idxn))
    PeriodicInterval(Interval(T(l),T(r)), support(dict))
end

period(b::WaveletBasis) = domaintype(b)(1)

interpolation_grid(b::WaveletBasis{T}) where {T} = DyadicPeriodicEquispacedGrid(length(b), support(b))


ordering(b::WaveletBasis{T,S,Wvl}) where {T,S<:Side} = wavelet_indices(dyadic_length(b))
ordering(b::WaveletBasis{T,S,Scl}) where {T,S<:Side} = scaling_indices(dyadic_length(b))


native_index(b::WaveletBasis{T,S,Wvl}, idx::Int) where {T,S<:Side} =
    wavelet_index(dyadic_length(b), idx)
linear_index(b::WaveletBasis{T,S,Wvl}, idxn::WaveletIndex) where {T,S<:Side} =
    value(idxn)
native_index(b::WaveletBasis{T,S,Scl}, idx::Int) where {T,S<:Side} =
    scaling_index(dyadic_length(b), idx)
linear_index(b::WaveletBasis{T,S,Scl}, idxn::WaveletIndex) where {T,S<:Side} =
    scaling_value(idxn)

approx_length(::WaveletBasis, n) = 1<<round(Int, log2(size_l))

extensionsize(b::WaveletBasis) = 2*length(b)


function unsafe_eval_element(dict::WaveletBasis{T,S,K,scaled}, idxn::WaveletIndex, x; xtol=1e-4, options...) where {T,S,K,scaled}
    @warn "Wavelet evaluation in arbitrary point is slow for $(typeof(dict))" maxlog=3
    evaluate_periodic(S(), kind(idxn), wavelet(dict), level(idxn), offset(idxn), x; xtol = xtol, scaled=scaled, options...)
end

unsafe_eval_element1(dict::WaveletBasis, idxn::WaveletIndex, grid::DyadicPeriodicEquispacedGrid; options...) =
    _unsafe_eval_element_in_dyadic_grid(dict, idxn, grid; options...)

function unsafe_eval_element1(dict::WaveletBasis, idxn::WaveletIndex, grid::PeriodicEquispacedGrid; options...)
    if isdyadic(length(grid))
        _unsafe_eval_element_in_dyadic_grid(dict, idxn, grid; options...)
    else
        error("No algorithm available to evaluate $(typeof(dict)) in $(typeof(grid))")
        _default_unsafe_eval_element_in_grid(dict, idxn, grid)
    end
end

function _unsafe_eval_element_in_dyadic_grid(dict::WaveletBasis{T,S,K,scaled}, idxn::WaveletIndex, grid::AbstractGrid; options...) where {T,S,K,scaled}
    @assert isdyadic(length(grid))
    evaluate_periodic_in_dyadic_points(S(), kind(idxn), wavelet(dict), level(idxn), offset(idxn), Int(log2(length(grid)));scaled=scaled)
end
DiscreteWaveletTransform(dict::WaveletBasis{T,S,Scl,scaled}) where {T,S,scaled} =
    DiscreteWaveletTransform(waveletbasis(dict), dict, wavelet(dict), side(dict), dyadic_length(dict), scaled)

DiscreteWaveletTransform(dict::WaveletBasis{T,S,Wvl,scaled}) where {T,S,scaled} =
    DiscreteWaveletTransform(dict, scalingbasis(dict), wavelet(dict), side(dict), dyadic_length(dict), scaled)

InverseDiscreteWaveletTransform(dict::WaveletBasis{T,S,Scl,scaled}) where {T,S,scaled} =
    InverseDiscreteWaveletTransform(dict, waveletbasis(dict), wavelet(dict), side(dict), dyadic_length(dict), scaled)

InverseDiscreteWaveletTransform(dict::WaveletBasis{T,S,Wvl,scaled}) where {T,S,scaled} =
    InverseDiscreteWaveletTransform(scalingbasis(dict), dict, wavelet(dict), side(dict), dyadic_length(dict), scaled)

DiscreteWaveletTransform(dict::WaveletTensorProductDict) =
    TensorProductOperator(map(DiscreteWaveletTransform, elements(dict))...)
InverseDiscreteWaveletTransform(dict::WaveletTensorProductDict) =
    TensorProductOperator(map(InverseDiscreteWaveletTransform, elements(dict))...)


evaluation(::Type{T}, s::WaveletBasis, dgs::GridBasis, grid::AbstractGrid; options...) where {T} =
    error("No algorithm available to evaluate $(typeof(s)) in $(typeof(grid))")

evaluation(::Type{T}, dict::WaveletBasis{T,S,Wvl}, dgs::GridBasis, grid::DyadicPeriodicEquispacedGrid; options...) where {T,S} =
    evaluation(T, scalingbasis(dict), dgs, grid; options...)*InverseDiscreteWaveletTransform(dict)

function evaluation(::Type{T}, dict::WaveletBasis{T,S,Scl,scaled}, dgs::GridBasis, grid::DyadicPeriodicEquispacedGrid; options...) where {T,S,scaled}
    d = dyadic_length(grid)
    w = wavelet(dict)
    j = dyadic_length(dict)
    s = side(dict)
    coefs = zeros(dict)
    coefs[1] = 1
    y = evaluate_periodic_scaling_basis_in_dyadic_points(s, w, coefs, d; scaled=scaled)
    x = dyadicpointsofrecursion(s, scaling, w, j, 0, d)
    offset = findfirst(x.==0)-1
    if length(x)>=length(y)
        return VerticalBandedOperator{T}(dict, dgs, y, 1<<(d-j), offset)
    end
    a = vcat(y[end-offset+1:end],y[1:length(x)-offset])
    VerticalBandedOperator{T}(dict, dgs, a, 1<<(d-j), offset)
end

evaluation(::Type{T}, s::WaveletBasis{T}, dgs::GridBasis, grid::PeriodicEquispacedGrid; options...) where {T} =
    isdyadic(grid) ?
        evaluation(T, s, dgs, DyadicPeriodicEquispacedGrid(grid); options...) :
        error("No algorithm available to evaluate $(typeof(s)) in $(typeof(grid))")

function evaluation_matrix!(a::AbstractMatrix, dict::WaveletBasis{T,S,K,scaled}, pts::DyadicPeriodicEquispacedGrid) where {T,S,K,scaled}
    @assert size(a,1) == length(pts)
    @assert size(a,2) == length(dict)

    s = side(dict)
    w = wavelet(dict)
    d = dyadic_length(pts)

    f = zeros(length(pts))
    SS = EvalPeriodicScratchSpace(s, w, dyadic_length(dict), d)
    for index in ordering(dict)
        evaluate_periodic_in_dyadic_points!(f, s, kind(index), w, level(index), offset(index), d, SS; scaled=scaled)
        for i in 1:length(f)
            a[i,value(index)] = f[i]
        end
    end
    a
end

export BiorthogonalWaveletBasis
abstract type BiorthogonalWaveletBasis{T,S,K,scaled} <: WaveletBasis{T,S,K,scaled} end

isbasis(b::BiorthogonalWaveletBasis) = true

export OrthogonalWaveletBasis
abstract type OrthogonalWaveletBasis{T,S,K,scaled} <: BiorthogonalWaveletBasis{T,S,K,scaled} end

isbasis(b::OrthogonalWaveletBasis) = true
isorthogonal(b::OrthogonalWaveletBasis, ::FourierMeasure) = true
isorthonormal(b::OrthogonalWaveletBasis, ::FourierMeasure) = true

export DaubechiesWaveletBasis
struct DaubechiesWaveletBasis{P,T,S,K,scaled} <: OrthogonalWaveletBasis{T,S,K,scaled}
    w   ::    DaubechiesWavelet{P,T}
    L   ::    Int
end

export scalingbasis
"""
    scalingbasis(w::DiscreteWavelet, L::Int, ::Type{S}=Prl)
"""
scalingbasis(w::DaubechiesWavelet{P,T}, L::Int, ::Type{S}=Prl, scaled=true) where {P,T,S} =
    Daubechiesscalingbasis(P, L, T, scaled)
waveletbasis(w::DaubechiesWavelet{P,T}, L::Int, ::Type{S}=Prl, scaled=true) where {P,T,S} =
    DaubechiesWaveletBasis(P, L, T, scaled)
scalingbasis(dict::DaubechiesWaveletBasis{P,T,S,K,scaled}) where {P,T,S,K,scaled} =
    DaubechiesWaveletBasis{P,T,S,Scl,scaled}(dict.w, dict.L)

export waveletbasis
"""
    waveletbasis(w::DiscreteWavelet, L::Int, ::Type{S}=Prl)
"""
waveletbasis(dict::DaubechiesWaveletBasis{P,T,S,K,scaled}) where {P,T,S,K,scaled} =
    DaubechiesWaveletBasis{P,T,S,Wvl,scaled}(dict.w, dict.L)

DaubechiesWaveletBasis(P::Int, L::Int, ::Type{T} = Float64, scaled::Bool=true) where {T<:Number} =
    DaubechiesWaveletBasis{P,T,Prl,Wvl,scaled}(DaubechiesWavelet{P,T}(), L)

Daubechiesscalingbasis(P::Int, L::Int, ::Type{T} = Float64, scaled::Bool=true) where {T<:Number} =
    DaubechiesWaveletBasis{P,T,Prl,Scl,scaled}(DaubechiesWavelet{P,T}(), L)

promote_domaintype(b::DaubechiesWaveletBasis{P,T,SIDE,KIND,scaled}, ::Type{S}) where {P,T,S,SIDE,KIND,scaled} =
    DaubechiesWaveletBasis{P,promote_type(T,S),SIDE,KIND,scaled}(DaubechiesWavelet{P,promote_type(T,S)}(), dyadic_length(b))

iscompatible(src1::DaubechiesWaveletBasis{P,T1,S1,K1,scaled}, src2::DaubechiesWaveletBasis{P,T2,S2,K2,scaled}) where {P,T1,T2,S1,S2,K1,K2,scaled} = true

# Note no check on existence of CDFXY is present.
export CDFWaveletBasis
struct CDFWaveletBasis{P,Q,T,S,K,scaled} <: BiorthogonalWaveletBasis{T,S,K,scaled}
    w   ::    CDFWavelet{P,Q,T}
    L   ::    Int
end

scalingbasis(w::CDFWavelet{P,Q,T}, L::Int, ::Type{S}=Prl, scaled::Bool=true) where {P,Q,T,S} =
    CDFscalingbasis(P, Q, L, S, T, scaled)
waveletbasis(w::CDFWavelet{P,Q,T}, L::Int, ::Type{S}=Prl, scaled::Bool=true) where {P,Q,T,S} =
    CDFWaveletBasis(P, Q, L, S, T, scaled)

CDFWaveletBasis(P::Int, Q::Int, L::Int, ::Type{S}=Prl, ::Type{T} = Float64, scaled::Bool=true) where {T,S<:Side} =
    CDFWaveletBasis{P,Q,T,S,Wvl,scaled}(CDFWavelet{P,Q,T}(),L)

CDFscalingbasis(P::Int, Q::Int, L::Int, ::Type{S}=Prl, ::Type{T} = Float64, scaled::Bool=true) where {T,S<:Side} =
    CDFWaveletBasis{P,Q,T,S,Scl,scaled}(CDFWavelet{P,Q,T}(),L)

scalingbasis(dict::CDFWaveletBasis{P,Q,T,S,K,scaled}) where {P,Q,T,S,K,scaled} =
    CDFWaveletBasis{P,Q,T,S,Scl,scaled}(dict.w, dict.L)

waveletbasis(dict::CDFWaveletBasis{P,Q,T,S,K,scaled}) where {P,Q,T,S,K,scaled} =
    CDFWaveletBasis{P,Q,T,S,Wvl,scaled}(dict.w, dict.L)

promote_domaintype(b::CDFWaveletBasis{P,Q,T,SIDE,KIND,scaled}, ::Type{S}) where {P,Q,T,S,SIDE,KIND,scaled} =
    CDFWaveletBasis{P,Q,promote_type(T,S),SIDE,KIND,scaled}(CDFWavelet{P,Q,promote_type(T,S)}(), dyadic_length(b))

iscompatible(src1::CDFWaveletBasis{P,Q,T1,S1,K1,scaled}, src2::CDFWaveletBasis{P,Q,T2,S2,K2,scaled}) where {P,Q,T1,T2,S1,S2,K1,K2,scaled} = true

GenericPeriodicEquispacedTranslates(dict::CDFWaveletBasis{P,Q,T,Prl,Scl,scaled}) where {P,Q,T,scaled} =
    BSplineTranslatesBasis{T,P-1,scaled}(length(dict))
#     GenericPeriodicEquispacedTranslates(interpolation_grid(dict), x->sqrt(T(length(dict)))*evaluate_BSpline(Val(P-1),
#         length(dict)*x, T), Interval(T(0),T(P)/length(dict)))
# GenericPeriodicEquispacedTranslates(dict::CDFWaveletBasis{P,Q,T,Prl,Scl,false}) where {P,Q,T} =
#     GenericPeriodicEquispacedTranslates(interpolation_grid(dict), x->evaluate_BSpline(Val(P-1),
#         length(dict)*x, T), Interval(T(0),T(P)/length(dict)))


for GRID in (:DyadicPeriodicEquispacedGrid, :AbstractIntervalGrid, :PeriodicEquispacedGrid)
    @eval begin
        evaluation(::Type{T}, dict::CDFWaveletBasis{P,Q,T,Prl,Scl}, dgs::GridBasis, grid::$GRID; options...) where {P,Q,T} =
            evaluation(T, GenericPeriodicEquispacedTranslates(dict), dgs, grid; options...)
        evaluation(::Type{T}, dict::CDFWaveletBasis{P,Q,T,Prl,Wvl}, dgs::GridBasis, grid::$GRID; options...) where {P,Q,T} =
            evaluation(T, scalingbasis(dict), dgs, grid; options...)*InverseDiscreteWaveletTransform(dict)

        unsafe_eval_element1(dict::CDFWaveletBasis{P,Q,T,Prl,Scl}, idxn::WaveletIndex, grid::$GRID; options...) where {P,Q,T} =
            unsafe_eval_element1(GenericPeriodicEquispacedTranslates(dict), value(idxn), grid; options...)
        function unsafe_eval_element1(dict::CDFWaveletBasis{P,Q,T,Prl,Wvl}, idxn::WaveletIndex, grid::$GRID; options...) where {P,Q,T}
            e = zeros(dict)
            e[value(idxn)] = 1
            evaluation(dict, GridBasis(grid), grid; options...)*e
        end
    end
end

function unsafe_eval_element(dict::CDFWaveletBasis{P,Q,T,S,K,scaled}, idxn::WaveletIndex, x; xtol=1e-4, options...) where {P,Q,T,S,K,scaled}
    evaluate_periodic(S(), kind(idxn), wavelet(dict), level(idxn), offset(idxn), x; xtol = xtol, scaled=scaled, options...)
end

@recipe function f(F::WaveletBasis; plot_complex = false, n=200)
    legend --> false
    grid = plotgrid(F,n)
    for i in eachindex(F)
        @series begin
            vals = F[i](grid)
            grid, postprocess(F[i],grid,vals)
        end
    end
    nothing
end

plotgrid(b::WaveletBasis, n) = DyadicPeriodicEquispacedGrid(1<<round(Int,log2(n)), support(b))
plotgrid(b::CDFWaveletBasis{P,Q,T,Prl}, n) where {P,Q,T} = PeriodicEquispacedGrid(round(Int,n/length(b))*length(b), support(b))


# export DWTSamplingOperator
# """
#     struct DWTSamplingOperator <: SamplingOperator
#
# A `DWTSamplingOperator` is an operator that maps a function to scaling coefficients.
# """
# struct DWTSamplingOperator <: SamplingOperator
#     sampler :: GridSampling
#     weight  :: DictionaryOperator
#     scratch :: Array
#
# 	# An inner constructor to enforce that the operators match
# 	function DWTSamplingOperator(sampler::GridSampling, weight::DictionaryOperator{ELT}) where {ELT}
#         @assert real(ELT) == eltype(eltype(grid(sampler)))
#         @assert size(weight, 2) == length(grid(sampler))
# 		new(sampler, weight, zeros(src(weight)))
#     end
# end
# using WaveletsEvaluation.DWT: quad_sf_weights
# convert(::Type{OP}, dwt::DWTSamplingOperator) where {OP<:DictionaryOperator} = dwt.weight
# promote_rule(::Type{OP}, ::DWTSamplingOperator) where{OP<:DictionaryOperator} = OP

# """
# A `WeightOperator` is an operator that maps function values to scaling coefficients.
# """
# function WeightOperator(basis::WaveletBasis, oversampling::Int=1, recursion::Int=0)
#     wav = wavelet(basis)
#     @assert coefficienttype(basis) == eltype(wav)
#     WeightOperator(wav, oversampling, dyadic_length(basis), recursion)
# end
#
# function WeightOperator(wav::DiscreteWavelet{T}, oversampling::Int, j::Int, d::Int) where {T}
#     @assert isdyadic.(oversampling)
#     w = quad_sf_weights(Dual, scaling, wav, oversampling*support_length(Dual, scaling, wav), d)
#     # rescaling of weights because the previous step assumed sum(w)==1
#     w .= w ./ sqrt(T(1<<j))
#     src_size = 1<<(d+j+oversampling>>1)
#     step = 1<<(d+oversampling>>1)
#     os = mod(step*InfiniteVectors.offset(filter(Dual, scaling, wav))-1, src_size)+1
#     try
#         HorizontalBandedOperator(GridBasis(DyadicPeriodicEquispacedGrid(1<<(d+j+oversampling>>1), UnitInterval{T}())), scalingbasis(wav, j) , w, step, os)
#     catch y
#         if isa(y, AssertionError)
#             error("Support of dual wavelet basis exceeds the width of the domain. Try more elements in your basis.")
#         else
#             rethrow(y)
#         end
#     end
# end
#
# function DWTSamplingOperator(dict::Dictionary, oversampling::Int=1, recursion::Int=0)
#     weight = WeightOperator(dict, oversampling, recursion)
#     sampler = GridSampling(GridBasis(dwt_oversampled_grid(dict, oversampling, recursion)))
#     DWTSamplingOperator(sampler, weight)
# end
#
# dwt_oversampled_grid(dict::Dictionary, oversampling::Int, recursion::Int) =
#     GridArrays.extend(interpolation_grid(dict), 1<<(recursion+oversampling>>1))
#
# src_space(op::DWTSamplingOperator) = src_space(op.sampler)
# dest(op::DWTSamplingOperator) = dest(op.weight)
#
# apply(op::DWTSamplingOperator, f) = apply!(zeros(dest(op)), op, f)
# apply!(result, op::DWTSamplingOperator, f) = apply!(op.weight, result, apply!(op.scratch, op.sampler, f))

##################
# Tensor methods
##################
# WeightOperator(dict::WaveletTensorProductDict, oversampling::Vector{Int}, recursion::Vector{Int}) =
#     TensorProductOperator([WeightOperator(di, osi, reci) for (di, osi, reci) in zip(elements(dict), oversampling, recursion)]...)
# _weight(dict::WaveletTensorProductDict, dyadic_os) =
#     TensorProductOperator([_weight(di, dyadic_os) for di in elements(dict)]...)
evaluation(::Type{T}, s::WaveletTensorProductDict, dgs::GridBasis, grid::AbstractSubGrid; options...) where {T} =
    restriction(GridBasis{coefficienttype(s)}(supergrid(grid)), GridBasis{coefficienttype(s)}(grid))*
        evaluation(T, s, GridBasis{coefficienttype(s)}(supergrid(grid)), supergrid(grid))
evaluation(::Type{T}, s::WaveletTensorProductDict, dgs::GridBasis, grid::ProductGrid; options...) where {T} =
    TensorProductOperator([evaluation(T, dict, GridBasis{coefficienttype(s)}(g), g) for (dict, g) in zip(elements(s), elements(grid))]...)

export wavelet_dual
"""
    wavelet_dual(dict::WaveletBasis)

    Create a similar wavelet basis, but replacing the wavelet with its biorthogonal dual
"""
wavelet_dual(w::DaubechiesWaveletBasis{P,T,S,K,scaled}) where {P,T,S,K,scaled} =
    DaubechiesWaveletBasis{P,T,inv(S),K,scaled}(DaubechiesWavelet{P,T}(),dyadic_length(w))
wavelet_dual(w::CDFWaveletBasis{P,Q,T,S,K,scaled}) where {P,Q,T,S,K,scaled} =
    CDFWaveletBasis{P,Q,T,inv(S),K,scaled}(CDFWavelet{P,Q,T}(),dyadic_length(w))

# function _weight(dest::WaveletBasis, dyadic_os)
#     j = dyadic_length(dest)
#     wav = wavelet(dest)
#     if dyadic_os == 0 # Just a choice I made.
#         return WeightOperator(wav, 1, j, 0)
#     else
#         return WeightOperator(wav, 2, j, dyadic_os-1)
#     end
# end
export LevelWeightingOperator
"""
    struct LevelWeightingOperator{T} <: DictionaryOperator{T}

An operator that weights the coefficients of the different wavelet levels differently.

# Example
```jldocs
julia> reshape(diag(LevelWeightingOperator(DaubechiesWaveletBasis(2,2)^2,[1.,.1])),4,4)
4×4 Array{Float64,2}:
 1.0  1.0  0.1  0.1
 1.0  1.0  0.1  0.1
 0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1
```
"""
struct LevelWeightingOperator{T} <: DictionaryOperator{T}
    src::WaveletTensorProductDict
    weights::Vector{T}
    function LevelWeightingOperator{T}(src::WaveletTensorProductDict,weights::Vector{T}) where T
        w = ones(T,minimum(map(dyadic_length,elements(src))))
        @assert minimum(map(dyadic_length,elements(src))) >= length(weights)
        copyto!(w,length(w)-length(weights)+1,weights)
        new{T}(src,w)
    end
    LevelWeightingOperator{T}(src::WaveletTensorProductDict,weights::Vector{S}) where {T,S} =
        LevelWeightingOperator{T}(src, convert.(T,weights))
    LevelWeightingOperator(src::WaveletTensorProductDict,weights::Vector{T}) where {T} =
        LevelWeightingOperator{T}(src, weights)
end

isinplace(::LevelWeightingOperator) = true
isdiag(::LevelWeightingOperator) = true
dest(op::LevelWeightingOperator) = src(op)
diag(op::LevelWeightingOperator) = reshape(op*ones(src(op)),length(src(op)))
function apply_inplace!(op::LevelWeightingOperator, src)
    limits = [i>>j for i in size(op.src),  j in length(op.weights):-1:1]
    limits[:,1] .+= 1
    for i in CartesianIndices(size(src))
        src[i] .*= op.weights[selectindex(i, limits)]
    end
    src
end

function selectindex(index::CartesianIndex, limits::Matrix{Int})
    r = 1
    for l in CartesianIndices(size(limits))
        if index[l.I[1]] > limits[l]
            r = l.I[2]
        end
    end
    r
end
adjoint(op::LevelWeightingOperator) = op

end

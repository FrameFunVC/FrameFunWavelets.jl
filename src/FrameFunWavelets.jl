module FrameFunWavelets
using WaveletsEvaluation, BasisFunctions, RecipesBase, DomainSets,
    CompactTranslatesDict, GridArrays, InfiniteVectors

using DomainSets: endpoints

using BasisFunctions: VerticalBandedOperator, SamplingOperator, HorizontalBandedOperator
using CompactTranslatesDict: PeriodicInterval

using WaveletsEvaluation.DWT: Side, Kind, Prl, Dul, Scl, Wvl, DiscreteWavelet, WaveletIndex
using WaveletsEvaluation.DWT: Filterbank, WaveletBoundary, SFilterBank
using WaveletsEvaluation.DWT: wavelet_index, scaling_indices, scaling_index
using WaveletsEvaluation.DWT: evaluate_periodic_in_dyadic_points!, dwt!, idwt!, offset,
    evaluate_periodic_scaling_basis_in_dyadic_points, _evaluate_periodic_scaling_basis_in_dyadic_points!


import WaveletsEvaluation.DWT: isdyadic, value, wavelet, kind

import GridArrays: similargrid, PeriodicEquispacedGrid

import BasisFunctions: subdict, name, interpolation_grid
import BasisFunctions: length, support, native_index, checkbounds, apply!
import BasisFunctions: hasextension, approx_length, resize, period, grid
import BasisFunctions: ordering, linear_index, hastransform, extension_size
import BasisFunctions: approximate_native_size
import BasisFunctions: transform_from_grid, transform_to_grid, unsafe_eval_element, unsafe_eval_element1
import BasisFunctions: isbasis, isorthogonal, evaluation_matrix!, gramoperator
import BasisFunctions: promote_domaintype, instantiate, iscompatible, plotgrid
import BasisFunctions: apply, dest, src_space


import Base:promote_eltype, ==, size

export DaubechiesWaveletBasis, CDFWaveletBasis, WaveletIndex, WaveletBasis, DaubechiesScalingBasis, CDFScalingBasis, ScalingBasis
export scaling_platform, WaveletTensorDict, dyadic_length, wavelet, kind, side, DyadicPeriodicEquispacedGrid


"""
A dyadic periodic equispaced grid is an equispaced grid that omits the right
endpoint and has length `n = 2^l`.
It has step `(b-a)/n`.
"""
struct DyadicPeriodicEquispacedGrid{T} <: AbstractEquispacedGrid{T}
    range   :: LinRange{T}
    l   ::  Int
    a   ::  T
    b   ::  T

    DyadicPeriodicEquispacedGrid{T}(n::Int, a, b) where {T} =
        (@assert isdyadic(n);new(LinRange(T(a),T(b),n+1)[1:end-1],Int(log2(n)),T(a),T(b)))
end


name(g::DyadicPeriodicEquispacedGrid) = "Dyadic periodic equispaced grid"
support(grid::DyadicPeriodicEquispacedGrid) = Interval(grid.a, grid.b)
isperiodic(::DyadicPeriodicEquispacedGrid) = true

dyadic_length(g::DyadicPeriodicEquispacedGrid) = g.l
length(g::DyadicPeriodicEquispacedGrid) = length(g.range)

DyadicPeriodicEquispacedGrid(n::Int, d::AbstractInterval) =
    DyadicPeriodicEquispacedGrid(n, endpoints(d)...)
similargrid(grid::DyadicPeriodicEquispacedGrid, ::Type{T}, n::Int) where {T} =
    DyadicPeriodicEquispacedGrid{T}(n, map(T, endpoints(support(grid)))...)
rescale(grid::DyadicPeriodicEquispacedGrid, a, b) =
    DyadicPeriodicEquispacedGrid{promote_type(typeof(a/2),typeof(b/2),eltype(grid))}(length(grid), a, b)
DyadicPeriodicEquispacedGrid(n::Int, a, b) =
    DyadicPeriodicEquispacedGrid{promote_type(typeof(a/2),typeof(b/2))}(n, a, b)
mapped_grid(grid::DyadicPeriodicEquispacedGrid, map::AffineMap) =
    DyadicPeriodicEquispacedGrid(length(grid), endpoints(map*support(grid))...)

isdyadic(::DyadicPeriodicEquispacedGrid) = true
isdyadic(grid::PeriodicEquispacedGrid) = isdyadic(length(grid))
isdyadic(::AbstractIntervalGrid) = false

PeriodicEquispacedGrid(g::DyadicPeriodicEquispacedGrid{T}) where {T} = PeriodicEquispacedGrid{T}(length(g), g.a, g.b)
DyadicPeriodicEquispacedGrid(g::PeriodicEquispacedGrid{T}) where {T} =
    (@assert isdyadic(g);DyadicPeriodicEquispacedGrid(length(g),g.a,g.b))

_extension_size(::DyadicPeriodicEquispacedGrid, n::Int, factor::Int) = (@assert isdyadic(factor);factor*n)
hasextension(::DyadicPeriodicEquispacedGrid) = true
extend(grid::DyadicPeriodicEquispacedGrid, factor::Int) =
    resize(grid, _extension_size(grid, length(grid), factor))

similargrid(g::DyadicPeriodicEquispacedGrid, a, b, T) = DyadicPeriodicEquispacedGrid{T}(g.l, a, b)
resize(g::DyadicPeriodicEquispacedGrid, n::Int) = DyadicPeriodicEquispacedGrid(n, g.a, g.b)


# We need this basic definition, otherwise equality does not seem to hold when T is BigFloat...
==(g1::DyadicPeriodicEquispacedGrid, g2::DyadicPeriodicEquispacedGrid) =
    (g1.l == g2.l) && (g1.a == g2.a) && (g1.b==g2.b)
==(g1::DyadicPeriodicEquispacedGrid, g2::PeriodicEquispacedGrid) =
    (length(g1) == length(g2)) && (g1.a == g2.a) && (g1.b==g2.b)

abstract type WaveletBasis{T,S,K} <: Dictionary1d{T,T} where {S <: Side,K<:Kind}
end

native_index(dict::WaveletBasis, idx::WaveletIndex) = idx


checkbounds(::Type{Bool}, dict::WaveletBasis, i::WaveletIndex) =
    checkbounds(Bool, dict, linear_index(dict, i))
"""
The number of levels in the wavelet basis
"""
dyadic_length(dict::WaveletBasis) = dict.L

length(dict::WaveletBasis) = 1<<dyadic_length(dict)
size(dict::WaveletBasis) = (length(dict),)

"""
The wavelet type
"""
WaveletsEvaluation.wavelet(dict::WaveletBasis) = dict.w

name(b::WaveletBasis) = "Basis of "*WaveletsEvaluation.DWT.name(wavelet(b))*" wavelets"

side(::WaveletBasis{T,S}) where {T,S} = S()

kind(::WaveletBasis{T,S,K}) where {T,S,K} = K()

# If only the first 2^L basis elements remain, this is equivalent to a smaller wavelet basis
function subdict(b::WaveletBasis, idx::OrdinalRange)
    if (step(idx)==1) && (first(idx) == 1) && isdyadic(last(idx))
        resize(b, last(idx))
    else
        LargeSubdict(b,idx)
    end
end

hasextension(b::WaveletBasis{T,S}) where{T,S} = true

approx_length(b::WaveletBasis, n::Int) = 1<<ceil(Int, log2(n))

resize(b::B, n::Int) where {B<:WaveletBasis} = B(wavelet(b),round(Int, log2(n)))

hasinterpolation_grid(::WaveletBasis) = true

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

support_length_of_compact_function(s::WaveletBasis{T,S,Scl}) where {T,S,Scl} = T(support_length(side(s), kind(s), wavelet(s)))/T(length(s))

approximate_native_size(::WaveletBasis, size_l) = 1<<ceil(Int, log2(size_l))

approx_length(::WaveletBasis, n) = 1<<round(Int, log2(size_l))

extension_size(b::WaveletBasis) = 2*length(b)


function unsafe_eval_element(dict::WaveletBasis{T,S}, idxn::WaveletIndex, x; xtol=1e-4, options...) where {T,S}
    evaluate_periodic(S(), kind(idxn), wavelet(dict), level(idxn), offset(idxn), x; xtol = xtol, options...)
end

unsafe_eval_element1(dict::WaveletBasis, idxn::WaveletIndex, grid::DyadicPeriodicEquispacedGrid; options...) =
    _unsafe_eval_element_in_dyadic_grid(dict, idxn, grid; options...)

function unsafe_eval_element1(dict::WaveletBasis, idxn::WaveletIndex, grid::PeriodicEquispacedGrid; options...)
    if isdyadic(length(grid))
        _unsafe_eval_element_in_dyadic_grid(dict, idxn, grid; options...)
    else
        _default_unsafe_eval_element_in_grid(dict, idxn, grid)
    end
end

_unsafe_eval_element_in_dyadic_grid(dict::WaveletBasis{T,S}, idxn::WaveletIndex, grid::AbstractGrid; options...) where {T,S} =
    evaluate_periodic_in_dyadic_points(S(), kind(idxn), wavelet(dict), level(idxn), offset(idxn), round(Int,log2(length(grid))))


function transform_from_grid(src, dest::WaveletBasis{T,S,Scl}, grid; options...) where {T,S}
    @assert iscompatible(dest, grid)
    scaling_transfrom_from_grid(src, dest, grid; options...)
end

function transform_from_grid(src, dest::WaveletBasis{T,S,Wvl}, grid; options...) where {T,S}
    @assert iscompatible(dest, grid)
    DiscreteWaveletTransform(dest)*scaling_transfrom_from_grid(src, dest, grid; options...)
end

function scaling_transfrom_from_grid(src, dest::WaveletBasis, grid; options...)
    j = dyadic_length(dest)
    l = dyadic_length(grid)
    dyadic_os = l-j
    _weight_operator(dest, dyadic_os)
end


function _weight_operator(dest::WaveletBasis, dyadic_os)
    j = dyadic_length(dest)
    wav = wavelet(dest)
    if dyadic_os == 0 # Just a choice I made.
        return WeightOperator(wav, 1, j, 0)
    else
        return WeightOperator(wav, 2, j, dyadic_os-1)
    end
end

function transform_to_grid(src::WaveletBasis, dest, grid; options...)
    @assert iscompatible(src, grid)
    EvalOperator(src, dest, dyadic_length(grid); options...)
end

"""
Transformation of wavelet coefficients to function values.
"""
struct DWTEvalOperator{T} <: DictionaryOperator{T}
    src::WaveletBasis
    dest::GridBasis

    s::Side
    w::DiscreteWavelet
    fb::Filterbank
    j::Int
    d::Int
    bnd::WaveletBoundary

    f::Vector{T}
    f_scaled::Vector{T}
    coefscopy::Vector{T}
    coefscopy2::Vector{T}
end

function EvalOperator(dict::WaveletBasis{T,S,Wvl}, dgs::GridBasis, d::Int; options...) where {T,S}
    w = wavelet(dict)
    j = dyadic_length(dict)
    s = side(dict)
    fb = SFilterBank(s, w)

    # DWT.evaluate_in_dyadic_points!(f, s, scaling, w, j, 0, d, scratch)
    f = evaluate_in_dyadic_points(s, scaling, w, j, 0, d)
    f_scaled = similar(f)
    coefscopy = zeros(1<<j)
    coefscopy2 = similar(coefscopy)

    DWTEvalOperator{T}(dict, dgs, s, w, fb, j, d, perbound, f, f_scaled, coefscopy, coefscopy2)
end

function apply!(op::DWTEvalOperator, y, coefs; options...)
    idwt!(op.coefscopy, coefs, op.fb, op.bnd, op.j, op.coefscopy2)
    _evaluate_periodic_scaling_basis_in_dyadic_points!(y, op.f, op.s, op.w, op.coefscopy, op.j, op.d, op.f_scaled)
    y
end

function EvalOperator(dict::WaveletBasis{T,S,Scl}, dgs::GridBasis, d::Int; options...) where {T,S}
    w = wavelet(dict)
    j = dyadic_length(dict)
    s = side(dict)
    coefs = zeros(dict)
    coefs[1] = 1
    y = evaluate_periodic_scaling_basis_in_dyadic_points(s, w, coefs, d)
    a, offset = _get_array_offset(y)
    VerticalBandedOperator(dict, dgs, a, 1<<(d-j), offset-1)
end

function _get_array_offset(a)
    b = a.!=0
    f = findfirst(b)
    if f==1
        if b[end]
            f = findlast(.!b)+1
            L = sum(b)
            vcat(a[f:end],a[1:L-length(a)+f]), f
        else
            a[f:f+sum(b)-1], f
        end
    else
        a[f:f+sum(b)-1], f
    end
end

struct DiscreteWaveletTransform{T} <: DictionaryOperator{T}
    src::WaveletBasis{T,S,Scl} where S
    dest::WaveletBasis{T,S,Wvl} where S

    fb::Filterbank
    j::Int

    scratch::Vector{T}
end

DiscreteWaveletTransform(dict::WaveletBasis{T,S,Scl}) where{T,S} =
    DiscreteWaveletTransform(dict, WaveletBasis(dict))

DiscreteWaveletTransform(dict::WaveletBasis{T,S,Wvl}) where{T,S} =
    DiscreteWaveletTransform(ScalingBasis(dict), dict)

function DiscreteWaveletTransform(src::WaveletBasis{T,S,Scl}, dest::WaveletBasis{T,S,Wvl}) where {T,S}
    w = wavelet(src)
    j = dyadic_length(src)
    s = side(src)
    fb = SFilterBank(s, w)
    scratch = zeros(T, 1<<j)
    DiscreteWaveletTransform{T}(src, dest, fb, j, scratch)
end

apply!(op::DiscreteWaveletTransform, dest, src; options...) =
    dwt!(dest, src, op.fb, perbound, op.j, op.scratch)

grid_evaluation_operator(s::WaveletBasis, dgs::GridBasis, grid::AbstractGrid; options...) =
    default_evaluation_operator(s, dgs; options...)

grid_evaluation_operator(s::WaveletBasis, dgs::GridBasis, grid::DyadicPeriodicEquispacedGrid; options...) =
    EvalOperator(s, dgs, dyadic_length(grid); options...)

function grid_evaluation_operator(s::WaveletBasis, dgs::GridBasis, subgrid::AbstractSubGrid; options...)
    # We make no attempt if the set has no associated grid
    if hasinterpolation_grid(s)
        # Is the associated grid of the same type as the supergrid at hand?
        if typeof(grid(s)) == typeof(supergrid(subgrid))
            # It is: we can use the evaluation operator of the supergrid
            super_dgs = gridbasis(s, supergrid(subgrid))
            E = evaluation_operator(s, super_dgs; options...)
            R = restriction_operator(super_dgs, dgs; options...)
            R*E
        else
            default_evaluation_operator(s, dgs; options...)
        end
    else
        default_evaluation_operator(s, dgs; options...)
    end
end

function evaluation_matrix!(a::AbstractMatrix, dict::WaveletBasis, pts::DyadicPeriodicEquispacedGrid)
    @assert size(a,1) == length(pts)
    @assert size(a,2) == length(dict)

    s = side(dict)
    w = wavelet(dict)
    d = dyadic_length(pts)

    f = zeros(length(pts))
    SS = EvalPeriodicScratchSpace(s, w, dyadic_length(dict), d)
    for index in ordering(dict)
        evaluate_periodic_in_dyadic_points!(f, s, kind(index), w, level(index), offset(index), d, SS)
        for i in 1:length(f)
            a[i,value(index)] = f[i]
        end
    end
    a
end

abstract type BiorthogonalWaveletBasis{T,S,K} <: WaveletBasis{T,S,K} end

isbasis(b::BiorthogonalWaveletBasis) = true

abstract type OrthogonalWaveletBasis{T,S,K} <: BiorthogonalWaveletBasis{T,S,K} end

isbasis(b::OrthogonalWaveletBasis) = true
isorthogonal(b::OrthogonalWaveletBasis, ::FourierMeasure) = true
isorthonormal(b::OrthogonalWaveletBasis, ::FourierMeasure) = true


struct DaubechiesWaveletBasis{P,T,S,K} <: OrthogonalWaveletBasis{T,S,K}
    w   ::    DaubechiesWavelet{P,T}
    L   ::    Int
end

ScalingBasis(w::DaubechiesWavelet{P,T}, L::Int, ::Type{S}=Prl) where {P,T,S} =
    DaubechiesScalingBasis(P, L, T)

ScalingBasis(dict::DaubechiesWaveletBasis{P,T,S,K}) where {P,T,S,K} =
    DaubechiesWaveletBasis{P,T,S,Scl}(dict.w, dict.L)

WaveletBasis(dict::DaubechiesWaveletBasis{P,T,S,K}) where {P,T,S,K} =
    DaubechiesWaveletBasis{P,T,S,Wvl}(dict.w, dict.L)

DaubechiesWaveletBasis(P::Int, L::Int, ::Type{T} = Float64) where {T} =
    DaubechiesWaveletBasis{P,T,Prl,Wvl}(DaubechiesWavelet{P,T}(), L)

DaubechiesScalingBasis(P::Int, L::Int, ::Type{T} = Float64) where {T} =
    DaubechiesWaveletBasis{P,T,Prl,Scl}(DaubechiesWavelet{P,T}(), L)

promote_domaintype(b::DaubechiesWaveletBasis{P,T,SIDE,KIND}, ::Type{S}) where {P,T,S,SIDE,KIND} =
    DaubechiesWaveletBasis{P,promote_type(T,S),SIDE,KIND}(DaubechiesWavelet{P,promote_type(T,S)}(), dyadic_length(b))

instantiate(::Type{DaubechiesWaveletBasis}, n, ::Type{T}) where {T} = DaubechiesWaveletBasis(3, approx_length(n), T)

iscompatible(src1::DaubechiesWaveletBasis{P,T1,S1,K1}, src2::DaubechiesWaveletBasis{P,T2,S2,K2}) where {P,T1,T2,S1,S2,K1,K2} = true

# Note no check on existence of CDFXY is present.
struct CDFWaveletBasis{P,Q,T,S,K} <: BiorthogonalWaveletBasis{T,S,K}
    w   ::    CDFWavelet{P,Q,T}
    L   ::    Int
end

ScalingBasis(w::CDFWavelet{P,Q,T}, L::Int, ::Type{S}=Prl) where {P,Q,T,S} =
    CDFScalingBasis(P, Q, L, S, T)

CDFWaveletBasis(P::Int, Q::Int, L::Int, ::Type{S}=Prl, ::Type{T} = Float64) where {T,S<:Side} =
    CDFWaveletBasis{P,Q,T,S,Wvl}(CDFWavelet{P,Q,T}(),L)

CDFScalingBasis(P::Int, Q::Int, L::Int, ::Type{S}=Prl, ::Type{T} = Float64) where {T,S<:Side} =
    CDFWaveletBasis{P,Q,T,S,Scl}(CDFWavelet{P,Q,T}(),L)

ScalingBasis(dict::CDFWaveletBasis{P,Q,T,S,K}) where {P,Q,T,S,K} =
    CDFWaveletBasis{P,Q,T,S,Scl}(dict.w, dict.L)

WaveletBasis(dict::CDFWaveletBasis{P,Q,T,S,K}) where {P,Q,T,S,K} =
    CDFWaveletBasis{P,Q,T,S,Wvl}(dict.w, dict.L)

promote_domaintype(b::CDFWaveletBasis{P,Q,T,SIDE,KIND}, ::Type{S}) where {P,Q,T,S,SIDE,KIND} =
    CDFWaveletBasis{P,Q,promote_type(T,S),SIDE,KIND}(CDFWavelet{P,Q,promote_type(T,S)}(), dyadic_length(b))

instantiate(::Type{CDFWaveletBasis}, n, ::Type{T}) where {T} = CDFWaveletBasis(2, 4, approx_length(n), T)

iscompatible(src1::CDFWaveletBasis{P,Q,T1,S1,K1}, src2::CDFWaveletBasis{P,Q,T2,S2,K2}) where {P,Q,T1,T2,S1,S2,K1,K2} = true


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

plotgrid(b::WaveletBasis, n) = DyadicPeriodicEquispacedGrid(round(Int,log2(n)), support(b))

"""
A `DWTSamplingOperator` is an operator that maps a function to scaling coefficients.
"""
struct DWTSamplingOperator <: SamplingOperator
    sampler :: GridSampling
    weight  :: DictionaryOperator
    scratch :: Array

	# An inner constructor to enforce that the operators match
	function DWTSamplingOperator(sampler::GridSampling, weight::DictionaryOperator{ELT}) where {ELT}
        @assert real(ELT) == eltype(eltype(grid(sampler)))
        @assert size(weight, 2) == length(grid(sampler))
		new(sampler, weight, zeros(src(weight)))
    end
end
using WaveletsEvaluation.DWT: quad_sf_weights
Base.convert(::Type{OP}, dwt::DWTSamplingOperator) where {OP<:DictionaryOperator} = dwt.weight
Base.promote_rule(::Type{OP}, ::DWTSamplingOperator) where{OP<:DictionaryOperator} = OP

"""
A `WeightOperator` is an operator that maps function values to scaling coefficients.
"""
function WeightOperator(basis::WaveletBasis, oversampling::Int=1, recursion::Int=0)
    wav = wavelet(basis)
    @assert coefficienttype(basis) == eltype(wav)
    WeightOperator(wav, oversampling, dyadic_length(basis), recursion)
end

function WeightOperator(wav::DiscreteWavelet{T}, oversampling::Int, j::Int, d::Int) where {T}
    @assert isdyadic.(oversampling)
    w = quad_sf_weights(Dual, scaling, wav, oversampling*support_length(Dual, scaling, wav), d)
    # rescaling of weights because the previous step assumed sum(w)==1
    w .= w ./ sqrt(T(1<<j))
    src_size = 1<<(d+j+oversampling>>1)
    step = 1<<(d+oversampling>>1)
    os = mod(step*InfiniteVectors.offset(filter(Dual, scaling, wav))-1, src_size)+1
    try
        HorizontalBandedOperator(GridBasis(DyadicPeriodicEquispacedGrid(1<<(d+j+oversampling>>1), UnitInterval{T}())), ScalingBasis(wav, j) , w, step, os)
    catch y
        if isa(y, AssertionError)
            error("Support of dual wavelet basis exceeds the width of the domain. Try more elements in your basis.")
        else
            rethrow(y)
        end
    end
end

function DWTSamplingOperator(dict::Dictionary, oversampling::Int=1, recursion::Int=0)
    weight = WeightOperator(dict, oversampling, recursion)
    sampler = GridSampling(GridBasis(dwt_oversampled_grid(dict, oversampling, recursion)))
    DWTSamplingOperator(sampler, weight)
end

dwt_oversampled_grid(dict::Dictionary, oversampling::Int, recursion::Int) =
    extend(interpolation_grid(dict), 1<<(recursion+oversampling>>1))

src_space(op::DWTSamplingOperator) = src_space(op.sampler)
dest(op::DWTSamplingOperator) = dest(op.weight)

apply(op::DWTSamplingOperator, f) = apply!(zeros(dest(op)), op, f)
apply!(result, op::DWTSamplingOperator, f) = apply!(op.weight, result, apply!(op.scratch, op.sampler, f))

##################
# Tensor methods
##################
const WaveletTensorDict2d = TensorProductDict{2,Tuple{B1,B2}} where {B1<:WaveletBasis,B2<:WaveletBasis}
const WaveletTensorDict3d = TensorProductDict{3,Tuple{B1,B2,B3}} where {B1<:WaveletBasis,B2<:WaveletBasis,B3<:WaveletBasis}
const WaveletTensorDict4d = TensorProductDict{4,Tuple{B1,B2,B3,B4}} where {B1<:WaveletBasis,B2<:WaveletBasis,B3<:WaveletBasis,B4<:WaveletBasis}
const WaveletTensorDict = Union{TensorProductDict{N,NTuple{N,B}} where {N,B<:WaveletBasis},WaveletTensorDict2d,WaveletTensorDict3d,WaveletTensorDict4d}

WeightOperator(dict::WaveletTensorDict, oversampling::Vector{Int}, recursion::Vector{Int}) =
    TensorProductOperator([WeightOperator(di, osi, reci) for (di, osi, reci) in zip(elements(dict), oversampling, recursion)]...)
_weight_operator(dict::WaveletTensorDict, dyadic_os) =
    TensorProductOperator([_weight_operator(di, dyadic_os) for di in elements(dict)]...)
grid_evaluation_operator(s::WaveletTensorDict, dgs::GridBasis, grid::AbstractSubGrid; options...) =
    restriction_operator(gridbasis(supergrid(grid), coefficienttype(s)), gridbasis(grid, coefficienttype(s)))*grid_evaluation_operator(s, gridbasis(supergrid(grid), coefficienttype(s)), supergrid(grid))
grid_evaluation_operator(s::WaveletTensorDict, dgs::GridBasis, grid::ProductGrid; options...) =
    TensorProductOperator([grid_evaluation_operator(dict, gridbasis(g, coefficienttype(s)), g) for (dict, g) in zip(elements(s), elements(grid))]...)

export wavelet_dual
"""
    wavelet_dual(dict::WaveletBasis)

    Create a similar wavelet basis, but replacing the wavelet with its biorthogonal dual
"""
wavelet_dual(w::OrthogonalWaveletBasis) = w
wavelet_dual(w::CDFWaveletBasis{P,Q,T,S,K}) where {P,Q,T,S,K} =
    CDFWaveletBasis{P,Q,T,inv(S),K}(CDFWavelet{P,Q,T}(),dyadic_length(w))

end
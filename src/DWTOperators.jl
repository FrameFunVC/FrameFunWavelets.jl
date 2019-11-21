module DWTOperators

using FrameFun.BasisFunctions

using WaveletsEvaluation.DWT: Filterbank, Side, DiscreteWavelet, scaling,
    evaluate_in_dyadic_points, _evaluate_periodic_scaling_basis_in_dyadic_points!,
    SFilterBank, idwt!, dwt!, perbound

import Base: size, inv, adjoint, transpose
import LinearAlgebra: mul!
import FrameFun.BasisFunctions: ArrayOperator, symbol, default_matrix, Matrix

"""
    struct DWTMatrix{T} <: AbstractMatrix{T}

Matrix that represents the periodic Discrete Wavelet Transform.
"""
struct DWTMatrix{T,W<:DiscreteWavelet,S<:Side} <: AbstractMatrix{T}
    fb  ::  Filterbank
    j   ::  Int
    scratch::Vector{T}

    function DWTMatrix(wavelet::DiscreteWavelet{T}, side::Side, j::Int) where T
        fb = SFilterBank(side, wavelet)
        scratch = zeros(T, 1<<j)
        new{T,typeof(wavelet),typeof(side)}(fb, j, scratch)
    end
end

size(A::DWTMatrix) = (1<<A.j, 1<<A.j)
function mul!(dest::AbstractVector{T}, op::DWTMatrix{T}, src::AbstractVector{T}) where T
    dwt!(dest, src, op.fb, perbound, op.j, op.scratch)
end

"""
    struct iDWTMatrix{T} <: AbstractMatrix{T}

Matrix that represents the periodic Discrete Wavelet Transform.
"""
struct iDWTMatrix{T,W<:DiscreteWavelet,S<:Side} <: AbstractMatrix{T}
    fb::Filterbank
    j::Int
    scratch::Vector{T}

    function iDWTMatrix(wavelet::DiscreteWavelet{T}, side::Side, j::Int) where T
        fb = SFilterBank(side, wavelet)
        scratch = zeros(T, 1<<j)
        new{T,typeof(wavelet),typeof(side)}(fb, j, scratch)
    end
end

size(A::iDWTMatrix) = (1<<A.j, 1<<A.j)
function mul!(dest::AbstractVector{T}, op::iDWTMatrix{T}, src::AbstractVector{T}) where T
    idwt!(dest, src, op.fb, perbound, op.j, op.scratch)
end


for op in (:inv, :adjoint, :transpose)
    @eval $op(A::DWTMatrix{T,W,S}) where {T,W,S} = iDWTMatrix(W(), S(), A.j)
    @eval $op(A::iDWTMatrix{T,W,S}) where {T,W,S} = DWTMatrix(W(), S(), A.j)
end

"""
    struct DWTEvalMatrix{T} <: AbstractMatrix{T}

Matrix that represents the evaluation matrix of a periodic wavelet
"""
struct DWTEvalMatrix{T,W<:DiscreteWavelet,S<:Side} <: AbstractMatrix{T}
    idwt::iDWTMatrix{T,W,S}
    s::S
    w::W
    j::Int
    d::Int

    f::Vector{T}
    f_scaled::Vector{T}
    coefscopy::Vector{T}

    function DWTEvalMatrix(wavelet::DiscreteWavelet{T}, side::Side, j::Int, d::Int) where T
        f = evaluate_in_dyadic_points(side, scaling, wavelet, j, 0, d)
        f_scaled = similar(f)
        coefscopy = zeros(T,1<<j)
        new{T,typeof(wavelet),typeof(side)}(iDWTMatrix(wavelet, side, j), side, wavelet, j, d, f, f_scaled, coefscopy)
    end
end

size(A::DWTEvalMatrix) = (1<<A.d, 1<<A.j)
function mul!(dest::AbstractVector{T}, op::DWTEvalMatrix{T}, src::AbstractVector{T}) where T
    mul!(op.coefscopy, op.idwt, src)
    _evaluate_periodic_scaling_basis_in_dyadic_points!(dest, op.f, op.s, op.w, op.coefscopy, op.j, op.d, op.f_scaled)
    dest
end


export DWTEvalOperator
"""
    struct DWTEvalOperator{T} <: ArrayOperator{T}

Transformation of wavelet coefficients to function values.
"""
struct DWTEvalOperator{T,W,S} <: ArrayOperator{T}
    A   ::  DWTEvalMatrix{T,W,S}
    src ::  Dictionary
    dest::  Dictionary
end

function DWTEvalOperator(src::Dictionary, dest::Dictionary, wavelet::DiscreteWavelet, side::Side, j::Int, d::Int)
    @assert length(src) == (1<<j)
    @assert length(dest) == (1<<d)
    DWTEvalOperator(DWTEvalMatrix(wavelet, side, j, d), src, dest)
end

ArrayOperator(A::DWTEvalMatrix, src::Dictionary, dest::Dictionary) =
    DWTEvalOperator(A, src, dest)

symbol(op::DWTEvalOperator) = "DWT_A"

Matrix(op::DWTEvalOperator) = default_matrix(op)

export DiscreteWaveletTransform
"""
    struct DiscreteWaveletTransform{T} <: ArrayOperator{T}

Operator that transforms scaling coefficients to wavelet coefficients
"""
struct DiscreteWaveletTransform{T,W,S} <: ArrayOperator{T}
    A   ::  DWTMatrix{T,W,S}
    src ::  Dictionary
    dest::  Dictionary
end

function DiscreteWaveletTransform(src::Dictionary, dest::Dictionary, wavelet::DiscreteWavelet, side::Side, j::Int)
    @assert length(src) == length(dest) == (1<<j)
    DiscreteWaveletTransform(DWTMatrix(wavelet, side, j), src, dest)
end

ArrayOperator(A::DWTMatrix, src::Dictionary, dest::Dictionary) =
    DiscreteWaveletTransform(A, src, dest)

symbol(op::DiscreteWaveletTransform) = "DWT"

Matrix(op::DiscreteWaveletTransform) = default_matrix(op)

export InverseDiscreteWaveletTransform
"""
    struct InverseDiscreteWaveletTransform{T} <: ArrayOperator{T}

Operator that transforms wavelet coefficients to scaling coefficients
"""
struct InverseDiscreteWaveletTransform{T,W,S} <: ArrayOperator{T}
    A   ::  iDWTMatrix{T,W,S}
    src ::  Dictionary
    dest::  Dictionary
end

function InverseDiscreteWaveletTransform(src::Dictionary, dest::Dictionary, wavelet::DiscreteWavelet, side::Side, j::Int)
    @assert length(src) == length(dest) == (1<<j)
    InverseDiscreteWaveletTransform(iDWTMatrix(wavelet, side, j), src, dest)
end

ArrayOperator(A::iDWTMatrix, src::Dictionary, dest::Dictionary) =
    InverseDiscreteWaveletTransform(A, src, dest)

symbol(op::InverseDiscreteWaveletTransform) = "iDWT"

Matrix(op::InverseDiscreteWaveletTransform) = default_matrix(op)


end

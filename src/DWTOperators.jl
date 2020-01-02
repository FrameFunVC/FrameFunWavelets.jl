module DWTOperators


using FrameFun.BasisFunctions, SparseArrays

using WaveletsEvaluation.DWT: Filterbank, Side, DiscreteWavelet, scaling,
    evaluate_in_dyadic_points, _evaluate_periodic_scaling_basis_in_dyadic_points!,
    InverseSFilterBank, SFilterBank, idwt!, dwt!, perbound

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
    scaled::Bool

    function DWTMatrix(wavelet::DiscreteWavelet{T}, side::Side, j::Int, scaled=true) where T
        fb = SFilterBank(side, wavelet, scaled)
        scratch = zeros(T, 1<<j)
        new{T,typeof(wavelet),typeof(side)}(fb, j, scratch, scaled)
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
    scaled::Bool

    function iDWTMatrix(wavelet::DiscreteWavelet{T}, side::Side, j::Int, scaled=true) where T
        fb = InverseSFilterBank(side, wavelet, scaled)
        scratch = zeros(T, 1<<j)
        new{T,typeof(wavelet),typeof(side)}(fb, j, scratch, scaled)
    end
end

size(A::iDWTMatrix) = (1<<A.j, 1<<A.j)
function mul!(dest::AbstractVector{T}, op::iDWTMatrix{T}, src::AbstractVector{T}) where T
    idwt!(dest, src, op.fb, perbound, op.j, op.scratch)
end

inv(A::DWTMatrix{T,W,S}) where {T,W,S} = iDWTMatrix(W(), S(), A.j, A.scaled)
inv(A::iDWTMatrix{T,W,S}) where {T,W,S} = DWTMatrix(W(), S(), A.j, A.scaled)

for op in (:transpose,:adjoint)
    @eval $op(A::DWTMatrix{T,W,S}) where {T,W,S} = iDWTMatrix(W(), inv(S()), A.j, A.scaled)
    @eval $op(A::iDWTMatrix{T,W,S}) where {T,W,S} = DWTMatrix(W(), inv(S()), A.j, A.scaled)
end


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

function DiscreteWaveletTransform(src::Dictionary, dest::Dictionary, wavelet::DiscreteWavelet, side::Side, j::Int, scaled=true)
    @assert length(src) == length(dest) == (1<<j)
    DiscreteWaveletTransform(DWTMatrix(wavelet, side, j, scaled), src, dest)
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

function InverseDiscreteWaveletTransform(src::Dictionary, dest::Dictionary, wavelet::DiscreteWavelet, side::Side, j::Int, scaled=true)
    @assert length(src) == length(dest) == (1<<j)
    InverseDiscreteWaveletTransform(iDWTMatrix(wavelet, side, j, scaled), src, dest)
end

ArrayOperator(A::iDWTMatrix, src::Dictionary, dest::Dictionary) =
    InverseDiscreteWaveletTransform(A, src, dest)

symbol(op::InverseDiscreteWaveletTransform) = "iDWT"

Matrix(op::InverseDiscreteWaveletTransform) = default_matrix(op)

end


module sparseiDWTEs

using WaveletsEvaluation, InfiniteVectors, FrameFunWavelets, BasisFunctions, StaticArrays, SparseArrays

import SparseArrays: sparse

export sparseiDWTEmatrix
function sparseiDWTEmatrix(basis::Union{WaveletBasis,WaveletTensorProductDict}, ix::AbstractVector, nz_tol=0)
    matrix = sparseiDWTMatrix(basis)
    sparseiDWTE(matrix, ix, nz_tol)
end

sparseiDWTMatrix(basis::WaveletBasis) = SparseiDWTMatrix(basis)
sparseiDWTMatrix(basis::WaveletTensorProductDict) = SparseiDWTTensor(basis)

function iDWTfilters(basis::WaveletBasis)
    primal_pair = InverseDiscreteWaveletTransform(basis).A.fb.primal_pair
    filters(primal_pair.f2, primal_pair.f1, dyadic_length(basis))
end

import Base: *, getindex, size
*(a::CompactPeriodicInfiniteVector) = a
function filters(g, h, J::Int)
    hs = [upsample(CompactPeriodicInfiniteVector(h,1<<(J-j)),1<<(j)) for j in (J-1):-1:0]
    gs = [upsample(CompactPeriodicInfiniteVector(g,1<<(J-j)),1<<(j)) for j in 0:J-1]

    r = Vector{CompactPeriodicInfiniteVector}(undef, J+1)
    r[1] = gs[1]
    t = hs[end]
    r[2] = gs[2]*t
    for j in 3:J
        t = hs[end-j+2]*t
        r[j] = gs[j]*t
    end
    r[end] = hs[1]*t
    r
end
struct SparseiDWTMatrix{T} <: AbstractMatrix{T}
    filters::Vector{CompactPeriodicInfiniteVector{T}}
    J::Int
    SparseiDWTMatrix(basis::WaveletBasis) =
        new{coefficienttype(basis)}(iDWTfilters(basis),dyadic_length(basis))
end
size(A::SparseiDWTMatrix) = (1<<A.J,1<<A.J)
function getindex(A::SparseiDWTMatrix, j::Int, i::Int)
    if i==1
        return A.filters[end][j]
    end
    band, offset = WaveletsEvaluation.DWT.wavelet_wavelet_index(1<<A.J, i)
    step = 1<<(A.J-band)
    A.filters[end-band-1][-step*offset+j-1]
end
using GridArrays.ModCartesianIndicesBase: ModUnitRange
function rowsupport!(vec::AbstractVector{Int}, M::SparseiDWTMatrix, row::Int)
    vec[1] = 1
    vec[2] = 2
    bandwidth = cld(InfiniteVectors.sublength(M.filters[1]),2)

    i = 3
    for band in 1:(M.J-1)
        step = 1<<(M.J-band)
        nzcolsinband = -floor(Int,(last(eachnonzeroindex(M.filters[end-band-1]))+1-row)/step):-ceil(Int,(first(eachnonzeroindex(M.filters[end-band-1]))  +1-row)/step)
        mod_nzcolsinband = ModUnitRange( 1<<band,nzcolsinband .+ 1)
        for coli in  mod_nzcolsinband
            vec[i] =  1<<band + coli
            i += 1
        end
    end
    i-1
end
function bandwidth(M::SparseiDWTMatrix)
    j = 2
    for i in 1:(M.J-1)
        j += min(cld(InfiniteVectors.sublength(M.filters[i]),1<<(i)),1<<(M.J-i))
    end
     j
end

struct SparseiDWTTensor{T,N} <: AbstractMatrix{T}
    SparseiDWTMatrices::SVector{N,SparseiDWTMatrix{T}}

end

SparseiDWTTensor(basis::WaveletTensorProductDict) =
    SparseiDWTTensor{coefficienttype(basis),dimension(basis)}(map(SparseiDWTMatrix,elements(basis)))

function size(A::SparseiDWTTensor)
    sizes = map(size, A.SparseiDWTMatrices)
    (prod(map(x->x[1],sizes)),prod(map(x->x[2],sizes)))
end

function getindex(A::SparseiDWTTensor, i::Int, j::Int)
    I = CartesianIndices(tuple(map(x->size(x,1),A.SparseiDWTMatrices)...))[i]
    J = CartesianIndices(tuple(map(x->size(x,2),A.SparseiDWTMatrices)...))[j]
    prod(map(getindex, A.SparseiDWTMatrices, I.I, J.I))
end
# TODO: do something more smart
sparseiDWTE(M::SparseiDWTMatrix, indices::AbstractVector, nz_tol=0) =
        sparseiDWTE(SparseiDWTTensor{eltype(M),1}(SVector{1,typeof(M)}(tuple(M))), map(CartesianIndex, indices), nz_tol)

function sparseiDWTE(M::SparseiDWTTensor{T,N}, indices::AbstractVector{CartesianIndex{N}}, nz_tol=0) where {N,T}
    b_supportlengths = map(bandwidth, M.SparseiDWTMatrices)
    b_supportlength = prod(b_supportlengths)
    nnz = b_supportlength*length(indices)

    colptr = Vector{Int}(undef, length(indices)+1)
    nzvals = Vector{T}(undef, nnz)
    rowvals = Vector{Int}(undef, nnz)

    rowvalscol = Vector{Int}(undef, b_supportlength)
    nzvalscol = Vector{T}(undef, b_supportlength)
    rowsupportcol = Matrix{Int}(undef, b_supportlength, N)
    colix = Vector{Int}(undef, b_supportlength)
    rowsupport_lengths = Vector{Int}(undef, N)
    colptr[1] = 1
    nzvalindex = 1
    L = LinearIndices(ntuple(k->size(M.SparseiDWTMatrices[k],1),Val(N)))
    for (i,k) in enumerate(indices)
        # support of element with index k
        for d in 1:N
            rowsupport_lengths[d] = rowsupport!(view(rowsupportcol,:,d), M.SparseiDWTMatrices[d], k.I[d])
        end
        colptrcol = 0
        row_support_indices = CartesianIndices(ntuple(d->rowsupport_lengths[d],Val(N)))
        row_support_length = length(row_support_indices)
        for j in row_support_indices
            l = CartesianIndex(ntuple(d->rowsupportcol[j[d],d], Val(N)))
            Mi = prod(ntuple(d->M.SparseiDWTMatrices[d][k[d],l[d]], Val(N)))
            if Mi != 0
                colptrcol += 1
                rowvalscol[colptrcol] = L[l]
                nzvalscol[colptrcol] = Mi
            end
        end
        for iii in colptrcol+1:b_supportlength
            rowvalscol[iii] = L[end]+1
        end

        for j in 1:b_supportlength
            colix[j] = j
        end
        sort!(colix ,1,row_support_length, InsertionSort,Base.Perm(Base.Order.ForwardOrdering(),rowvalscol))
        for j in 1:colptrcol
            rowvals[nzvalindex] = rowvalscol[colix[j]]
            nzvals[nzvalindex] = nzvalscol[colix[j]]
            nzvalindex += 1
        end
        colptr[i+1] = colptr[i] + colptrcol
    end
    nz_tol == 0 ?
        SparseMatrixCSC(size(M,1),length(indices),colptr,resize!(rowvals,nzvalindex-1),resize!(nzvals,nzvalindex-1))' :
        droptol!(SparseMatrixCSC(size(M,1),length(indices),colptr,resize!(rowvals,nzvalindex-1),resize!(nzvals,nzvalindex-1)),nz_tol)'
end



function sparse(A::InverseDiscreteWaveletTransform)
    matrix = sparseiDWTMatrix(dest(A))
    op = ArrayOperator(sparseiDWTE(matrix, 1:length(dest(A)), 0), src(A), dest(A))
    @assert op ≈ A
    op
end


end

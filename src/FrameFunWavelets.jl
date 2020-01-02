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



module SparseMult
using SparseArrays, LinearAlgebra

using SparseArrays: estimate_mulsize, prefer_sort
using LinearAlgebra: Adjoint

function myspmatmul(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    mA, nA = size(A)
    nB = size(B, 2)
    nA == size(B, 1) || throw(DimensionMismatch())

    rowvalA = rowvals(A); nzvalA = nonzeros(A)
    rowvalB = rowvals(B); nzvalB = nonzeros(B)
    nnzC = max(estimate_mulsize(mA, nnz(A), nA, nnz(B), nB) * 11 ÷ 10, mA)
    @show nnzC, nB
    colptrC = Vector{Ti}(undef, nB+1)
    rowvalC = Vector{Ti}(undef, nnzC)
    nzvalC = Vector{Tv}(undef, nnzC)
    nzpercol = nnzC ÷ max(nB, 1)


    @inbounds begin
        ip = 1
        xb = fill(false, mA)
        for i in 1:nB  # over all cols of C
            if ip + mA - 1 > nnzC
                @show "extra"
                nnzC += max(mA, nnzC>>2)
                resize!(rowvalC, nnzC)
                resize!(nzvalC, nnzC)
            end
            colptrC[i] = ip0 = ip
            k0 = ip - 1
            @show(length((nzrange(B, i))))
            @time for jp in nzrange(B, i) # over all nz elements of col i of B
                nzB = nzvalB[jp]
                j = rowvalB[jp] # j is row of jpth nz element of col i of B
                for kp in nzrange(A, j) # over all nz elements of col j of A
                    nzC = nzvalA[kp] * nzB
                    k = rowvalA[kp]
                    if xb[k] # if we did already handle C[k,i]
                        nzvalC[k+k0] += nzC
                    else
                        nzvalC[k+k0] = nzC
                        xb[k] = true
                        rowvalC[ip] = k
                        ip += 1
                    end
                end
            end
            if ip > ip0
                if prefer_sort(ip-k0, mA)
                    # in-place sort of indices. Effort: O(nnz*ln(nnz)).
                    sort!(rowvalC, ip0, ip-1, QuickSort, Base.Order.Forward)
                    for vp = ip0:ip-1
                        k = rowvalC[vp]
                        xb[k] = false
                        nzvalC[vp] = nzvalC[k+k0]
                    end
                else
                    # scan result vector (effort O(mA))
                    for k = 1:mA
                        if xb[k]
                            xb[k] = false
                            rowvalC[ip0] = k
                            nzvalC[ip0] = nzvalC[k+k0]
                            ip0 += 1
                        end
                    end
                end
            end
        colptrC[nB+1] = ip
    end
end
    resize!(rowvalC, ip - 1)
    resize!(nzvalC, ip - 1)

    # This modification of Gustavson algorithm has sorted row indices
    C = SparseMatrixCSC(mA, nB, colptrC, rowvalC, nzvalC)
    return C
end


function myspmatmul(A::SparseMatrixCSC{Tv,Ti}, B::Adjoint{Tv,SparseMatrixCSC{Tv,Ti}}) where {Tv,Ti}
    mA, nA = size(A)
    Bparent = B.parent
    nB = size(B, 2)
    nA == size(B, 1) || throw(DimensionMismatch())

    rowvalA = rowvals(A); nzvalA = nonzeros(A)
    rowvalBparent = rowvals(Bparent); nzvalBparent = nonzeros(Bparent)
    nnzC = max(estimate_mulsize(mA, nnz(A), nA, nnz(B), nB) * 11 ÷ 10, mA)

    colptrC = Vector{Ti}(undef, nB+1)
    rowvalC = Vector{Ti}(undef, nnzC)
    nzvalC = Vector{Tv}(undef, nnzC)
    nzpercol = nnzC ÷ max(nB, 1)


    @time begin
    @inbounds begin
        ip = 1
        xb = fill(false, nB) # all colomns of C
        for i in 1:nA  # over all cols of A
            if ip + nB - 1 > nnzC
                @show "extra"
                nnzC += max(mA, nnzC>>2)
                resize!(rowvalC, nnzC)
                resize!(nzvalC, nnzC)
            end
            # colptrC[i] = ip0 = ip
            # k0 = ip - 1

            for jp in nzrange(Bparent, i) # over all nz elements of row i of B

                nzB = nzvalBparent[jp]
                j = rowvalBparent[jp] # j is column of jpth nz element of row i of B

                for kp in nzrange(A, jp) # over all nz elements of col jp of A
                    nzC = nzvalA[kp]*nzB
                    k = rowvalA[kp]
                    if xb[k] # if we did already handle C[k,j]
                        nzvalC[k+k0] += nzC
                    else
                        nzvalC[k+k0] = nzC
                        xb[k] = true
                        rowvalC[ip] = k
                        ip += 1
                    end
                end
            end
            if ip > ip0
                if prefer_sort(ip-k0, mA)
                    # in-place sort of indices. Effort: O(nnz*ln(nnz)).
                    sort!(rowvalC, ip0, ip-1, QuickSort, Base.Order.Forward)
                    for vp = ip0:ip-1
                        k = rowvalC[vp]
                        xb[k] = false
                        nzvalC[vp] = nzvalC[k+k0]
                    end
                else
                    # scan result vector (effort O(mA))
                    for k = 1:mA
                        if xb[k]
                            xb[k] = false
                            rowvalC[ip0] = k
                            nzvalC[ip0] = nzvalC[k+k0]
                            ip0 += 1
                        end
                    end
                end
            end
        end
        colptrC[nB+1] = ip
    end
end
    resize!(rowvalC, ip - 1)
    resize!(nzvalC, ip - 1)

    # This modification of Gustavson algorithm has sorted row indices
    C = SparseMatrixCSC(mA, nB, colptrC, rowvalC, nzvalC)
    return C
end

end
end

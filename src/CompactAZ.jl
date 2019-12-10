module CompactAZ


module sparseiDWTEs

using WaveletsEvaluation, InfiniteVectors, FrameFunWavelets, BasisFunctions, StaticArrays, SparseArrays
export sparseiDWTEmatrix
function sparseiDWTEmatrix(basis::Union{WaveletBasis,WaveletTensorProductDict}, ix::AbstractVector)
    matrix = sparseiDWTMatrix(basis)
    sparseiDWTE(matrix, ix)
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

function sparseiDWTE(M::SparseiDWTMatrix{T}, indices::AbstractVector{Int}) where T
    b_supportlength = bandwidth(M)
    nnz = b_supportlength*length(indices)

    colptr = Vector{Int}(undef, length(indices)+1)
    nzvals = Vector{T}(undef, nnz)
    rowvals = Vector{Int}(undef, nnz)

    rowvalscol = Vector{Int}(undef, b_supportlength)
    nzvalscol = Vector{T}(undef, b_supportlength)
    rowsupportcol = Vector{Int}(undef, b_supportlength)
    colix = Vector{Int}(undef, b_supportlength)
    colptr[1] = 1
    nzvalindex = 1
    for (i,k) in enumerate(indices)
        # support of element with index k
        rowsupport_length = rowsupport!(rowsupportcol, M, k)
        colptrcol = 0
        for j in 1:rowsupport_length
            l = rowsupportcol[j]
            if M[k,l] != 0
                colptrcol += 1
                rowvalscol[colptrcol] = l
                nzvalscol[colptrcol] = M[k,l]
            end
        end

        for j in 1:rowsupport_length
            colix[j] = j
        end
        sort!(colix ,1,rowsupport_length, InsertionSort,Base.Perm(Base.Order.ForwardOrdering(),rowvalscol))
        for j in 1:colptrcol
            rowvals[nzvalindex] = rowvalscol[colix[j]]
            nzvals[nzvalindex] = nzvalscol[colix[j]]
            nzvalindex += 1
        end
        colptr[i+1] = colptr[i] + colptrcol
    end

    SparseMatrixCSC(1<<M.J,length(indices),colptr,resize!(rowvals,nzvalindex-1),resize!(nzvals,nzvalindex-1))'
end
struct SparseiDWTTensor{T,N} <: AbstractMatrix{T}
    SparseiDWTMatrices::SVector{N,SparseiDWTMatrix{T}}

    SparseiDWTTensor(basis::WaveletTensorProductDict) =
        new{coefficienttype(basis),dimension(basis)}(map(SparseiDWTMatrix,elements(basis)))
end
function size(A::SparseiDWTTensor)
    sizes = map(size, A.SparseiDWTMatrices)
    (prod(map(x->x[1],sizes)),prod(map(x->x[2],sizes)))
end

function getindex(A::SparseiDWTTensor, i::Int, j::Int)
    I = CartesianIndices(tuple(map(x->size(x,1),A.SparseiDWTMatrices)...))[i]
    J = CartesianIndices(tuple(map(x->size(x,2),A.SparseiDWTMatrices)...))[j]
    prod(map(getindex, A.SparseiDWTMatrices, I.I, J.I))
end

function sparseiDWTE(M::SparseiDWTTensor{T,N}, indices::AbstractVector{CartesianIndex{N}}) where {N,T}
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
            rowvalscol[iii] = L[end]
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

    SparseMatrixCSC(size(M,1),length(indices),colptr,resize!(rowvals,nzvalindex-1),resize!(nzvals,nzvalindex-1))'
end

end



using ...WaveletPlatforms, ...WaveletBases, FrameFun, BasisFunctions, FrameFunTranslates.CompactAZ.CompactFrameFunExtension, CompactTranslatesDict.CompactInfiniteVectors

using WaveletsEvaluation.DWT: Wvl,Scl
using FrameFunTranslates.CompactAZ.CompactFrameFunExtension: _nonzero_coefficients
import FrameFunTranslates.CompactAZ.CompactFrameFunExtension: ef_nonzero_coefficients, compactinfinitevectors, ef_sparse_reducedAAZAoperator, ef_sparseAZ_AAZAreductionsolver

function ef_nonzero_coefficients(samplingstyle::SamplingStyle, platform::ExtensionFramePlatform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, L; options...)
    os_grid = haskey(options, :os_grid) ? options[:os_grid] : oversampling_grid(samplingstyle, platform, param, L; options...)
    q = div.(L, 1 .<< param)
    spline_nonzero_coeffs = _nonzero_coefficients((compactsupport(samplingstyle, scalingplatform(platform.basisplatform), param,
        map(scalingplatform,platforms), supergrid(os_grid); options...)), q, mask(os_grid))
    op = InverseDiscreteWaveletTransform(dictionary(platform.basisplatform, param))'
    a = zeros(src(op))
    a[spline_nonzero_coeffs].=NaN
    findall(isnan.(op*a))
end

function compactinfinitevectors(ss::OversamplingStyle, bplatform::Platform, param, platforms::Tuple{<:AbstractWaveletPlatform{<:Number,Scl}}, os_grid::AbstractIntervalGrid; options...)
     tuple(compactinfinitevector(dictionary(bplatform, param), os_grid))
end
function compactinfinitevectors(ss::OversamplingStyle, bplatform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Scl}}}, os_grid::ProductGrid; options...)
     map(compactinfinitevector, map(dictionary, platforms, param), elements(os_grid))
end


function ef_sparseAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::ExtensionFramePlatform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, L, directsolver; verbose=false, options...)
    os_grid = haskey(options, :os_grid) ? options[:os_grid] : oversampling_grid(samplingstyle, platform, param, L; verbose=verbose, options...)
    rM = haskey(options, :sparse_reducedAAZAoperator) ? options[:sparse_reducedAAZAoperator] : sparse_reducedAAZAoperator(samplingstyle, platform, param, L; verbose=verbose, os_grid=os_grid, options...)
    verbose && @info "Sparse AZ: use $(directsolver) as solver for first sparse AZ step"
    FrameFunInterface.directsolver(rM; verbose=verbose, directsolver=directsolver, options...)
end


using FrameFunTranslates.CompactAZ.CompactFrameFunExtension: _sparseRAE, sparseidentity, overlappingindices, nonzero_rows, _nonzero_coefficients
using .sparseiDWTEs: sparseiDWTEmatrix
using SparseArrays
function ef_sparse_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::ExtensionFramePlatform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, L; verbose=false, nz_tol=0, options...)
    verbose && @info "SparseAZStyle: create sparse A-AZ^*A"
    os_grid = haskey(options, :os_grid) ? options[:os_grid] : oversampling_grid(samplingstyle, platform, param, L; options...)

    wavelet_frame = dictionary(platform, param)
    wavelet_basis = basis(wavelet_frame)

    scaling_platform = ExtensionFramePlatform(scalingplatform(platform.basisplatform), support(platform))
    scaling_platforms = scalingplatform.(platforms)
    scaling_param = 1 .<< param

    scaling_frame = dictionary(scaling_platform, param)
    scaling_basis = basis(scaling_frame)

    scaling_frame2 = azdual_dict(scaling_platform, param; samplingstyle=samplingstyle, L=L)
    scaling_basis2 = basis(scaling_frame2)

    cvecs = scaling_basis isa Dictionary1d ?
        tuple(compactinfinitevector(scaling_basis, supergrid(os_grid))) :
        map(compactinfinitevector, elements(scaling_basis), elements(supergrid(os_grid)))

    cvecs_dual =  scaling_basis2 isa Dictionary1d ?
        tuple(compactinfinitevector(scaling_basis2, supergrid(os_grid))) :
        map(compactinfinitevector, elements(scaling_basis2), elements(supergrid(os_grid)))

    scaling_nonzero_coefs = _nonzero_coefficients(compactsupport(cvecs), div.(L,scaling_param), mask(os_grid))
    ix1 = scaling_nonzero_coefs
    A = _sparseRAE(cvecs, os_grid, ix1, scaling_param; nz_tol=nz_tol, verbose=verbose, options...)
    nzrows = findall(nonzero_rows(A))
    ix2 = overlappingindices(cvecs_dual, os_grid, nzrows, scaling_param, L)
    ix3 = unique(sort(vcat(ix1,ix2)))
    # In principle we only need ix in next line, but we need the extension to ix3 somehow for the second next line
    Z = _sparseRAE(cvecs_dual, os_grid, ix3, scaling_param; nz_tol=nz_tol, verbose=verbose, options...)
    ImZA = sparseidentity(ix1,ix3)-Z'A

    RAE = _sparseRAE(cvecs, os_grid, ix3, scaling_param; nz_tol=nz_tol, verbose=verbose, options...)
    iDWTE = sparseiDWTEmatrix(wavelet_basis, scaling_nonzero_coefs)

    scaling_M = droptol!(RAE*ImZA,nz_tol)
    verbose && @info "Sparse AZ: removing everything smaller than $nz_tol"
    M = droptol!(scaling_M*iDWTE, nz_tol)
    verbose && @info "Sparse AZ: A-AZ^*A has size $(size(M)) and $(nnz(M)) nonzero elements ($(100nnz(M)/prod(size(M)))% fill)"

    # wavelet_nonzero_coefs = nonzero_rows(iDWTE')

    src = wavelet_frame#[wavelet_nonzero_coefs]
    dest = GridBasis{coefficienttype(src)}(os_grid)
    ArrayOperator(M, src, dest)
end
end

module CompactAZ

using FrameFun, BasisFunctions
using ...WaveletPlatforms, ...WaveletBases, FrameFunTranslates.CompactAZ.CompactFrameFunExtension

using WaveletsEvaluation.DWT: Wvl,Scl
using FrameFunTranslates.CompactAZ.CompactFrameFunExtension: _nonzero_coefficients, _nonzero_pointsindices, ef_true_nonzero_reducedAZ_AAZAreductionsolver, ef_true_nonzero_reducedAAZAoperator
import FrameFunTranslates.CompactAZ.CompactFrameFunExtension: ef_nonzero_coefficients, compactinfinitevectors, ef_nonzero_pointsindices,
    ef_sparse_reducedAAZAoperator, ef_sparseAZ_AAZAreductionsolver, ef_reducedAZ_AAZAreductionsolver, ef_reducedAAZAoperator, ef_AAZA_nonzero_row_indexset

using CompactTranslatesDict: compactinfinitevector

const WE_L = nonzero_coefficients
const WE_M = AAZA_nonzero_row_indexset

function ef_nonzero_coefficients(samplingstyle::SamplingStyle, platform::ExtensionFramePlatform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, L; options...)
    os_grid = haskey(options, :os_grid) ? options[:os_grid] : oversampling_grid(samplingstyle, platform, param, L; options...)
    q = div.(L, 1 .<< param)
    WE_K = _nonzero_coefficients((compactsupport(samplingstyle, scalingplatform(platform.basisplatform), param,
        map(scalingplatform,platforms), supergrid(os_grid); options...)), q, mask(os_grid))


    op = InverseDiscreteWaveletTransform(dictionary(platform.basisplatform, param))'
    a = zeros(src(op))
    a[WE_K].=NaN
    findall(isnan.(op*a))
end

function ef_nonzero_coefficients(samplingstyle::SamplingStyle, platform::ExtensionFramePlatform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Scl}}}, L; options...)
    os_grid = haskey(options, :os_grid) ? options[:os_grid] : oversampling_grid(samplingstyle, platform, param, L; options...)
    q = div.(L, 1 .<< param)
    WE_K = _nonzero_coefficients((compactsupport(samplingstyle, scalingplatform(platform.basisplatform), param,
        map(scalingplatform,platforms), supergrid(os_grid); options...)), q, mask(os_grid))
end

function ef_AAZA_nonzero_row_indexset(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, L; verbose=false, options...)
    os_grid = haskey(options, :os_grid) ? options[:os_grid] : oversampling_grid(samplingstyle, platform, param, L; options...)
    scaling_platform = ExtensionFramePlatform(scalingplatform(platform.basisplatform), support(platform))
    scaling_platforms = map(scalingplatform,platforms)
    verbose && @info "ReducedAZStyle: use nonzero cols of scaling basis to get nonzero rows"
    nonzero_coefs = ef_nonzero_coefficients(samplingstyle, scaling_platform, param, scaling_platforms, L; verbose=verbose,  options..., os_grid=os_grid)
    ef_AAZA_nonzero_row_indexset(samplingstyle, scaling_platform, param, scaling_platforms, L; options..., os_grid=os_grid, nonzero_coefs=nonzero_coefs)
end

function compactinfinitevectors(ss::OversamplingStyle, bplatform::Platform, param, platforms::Tuple{<:AbstractWaveletPlatform{<:Number,Scl}}, os_grid::AbstractIntervalGrid; options...)
     tuple(compactinfinitevector(dictionary(bplatform, param), os_grid))
end
function compactinfinitevectors(ss::OversamplingStyle, bplatform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Scl}}}, os_grid::ProductGrid; options...)
     map(compactinfinitevector, map(dictionary, platforms, param), elements(os_grid))
end

ef_nonzero_coefficients(ss::DiscreteStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Scl}}}, q, os_grid::AbstractGrid; options...) =
    _nonzero_coefficients(compactsupport(ss, platform.basisplatform, param, platforms, supergrid(os_grid); os_grid=supergrid(os_grid), options...), q, mask(os_grid))


ef_nonzero_pointsindices(ss::DiscreteStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Scl}}}, q, os_grid::AbstractGrid, ix, relative::Bool; options...) =
    _nonzero_pointsindices(compactsupport(ss, platform.basisplatform, param, platforms, supergrid(os_grid); os_grid=supergrid(os_grid), options...), q, mask(os_grid), ix, relative)

function ef_nonzero_pointsindices(ss::DiscreteStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, q, os_grid::AbstractGrid, ix, relative::Bool; options...)
    _nonzero_pointsindices(compactsupport(ss, scalingplatform(platform.basisplatform), param, map(scalingplatform,platforms), supergrid(os_grid); os_grid=supergrid(os_grid), options...), q, mask(os_grid), ix, relative)
end


ef_reducedAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform}}, L, directsolver; options...) =
    ef_true_nonzero_reducedAZ_AAZAreductionsolver(samplingstyle, platform, param, platforms, L, directsolver; options...)

ef_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::Platform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform}}, L; options...) =
    ef_true_nonzero_reducedAAZAoperator(samplingstyle, platform, param, platforms, L; options...)
using LinearAlgebra: Diagonal
function ef_sparseAZ_AAZAreductionsolver(samplingstyle::SamplingStyle, platform::ExtensionFramePlatform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, L, directsolver;
        verbose=false, weightedAZ=false, options...)
    os_grid = haskey(options, :os_grid) ? options[:os_grid] : oversampling_grid(samplingstyle, platform, param, L; verbose=verbose, options...)
    rM = haskey(options, :sparse_reducedAAZAoperator) ? options[:sparse_reducedAAZAoperator] : sparse_reducedAAZAoperator(samplingstyle, platform, param, L; verbose=verbose, os_grid=os_grid, options...)
    verbose && @info "Sparse AZ: use $(directsolver) as solver for first sparse AZ step"

    if weightedAZ
        verbose && @info "Weighted AZ"
        AZ_Cweight = haskey(options,:AZ_Cweight) ? options[:AZ_Cweight] : error("No options `AZ_Cweight`")
        AZ_Cweight*
            FrameFunInterface.directsolver(ArrayOperator(rM.A*Diagonal(diag(AZ_Cweight)),src(rM),dest(rM)); verbose=verbose, directsolver=directsolver, options...)
    else
        FrameFunInterface.directsolver(rM; verbose=verbose, directsolver=directsolver, options...)
    end
end


using FrameFunTranslates.CompactAZ.CompactFrameFunExtension: _sparseRAE, sparseidentity, overlappingindices, nonzero_rows, _nonzero_coefficients
using ...sparseiDWTEs: sparseiDWTEmatrix
using SparseArrays
function ef_sparse_reducedAAZAoperator(samplingstyle::SamplingStyle, platform::ExtensionFramePlatform, param, platforms::Tuple{Vararg{<:AbstractWaveletPlatform{<:Number,Wvl}}}, L; return_scalingandidwteoperators=false, verbose=false, nz_tol=0, options...)
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
    iDWTE = sparseiDWTEmatrix(wavelet_basis, scaling_nonzero_coefs, nz_tol)

    scaling_M = droptol!(RAE*ImZA,nz_tol)
    verbose && @info "Sparse AZ: removing everything smaller than $nz_tol"
    if return_scalingandidwteoperators
        return scaling_M, iDWTE
    end
    M = droptol!((scaling_M*iDWTE), nz_tol)
    verbose && @info "Sparse AZ: A-AZ^*A has size $(size(M)) and $(nnz(M)) nonzero elements ($(100nnz(M)/prod(size(M)))% fill)"

    # wavelet_nonzero_coefs = nonzero_rows(iDWTE')

    src = wavelet_frame#[wavelet_nonzero_coefs]
    dest = GridBasis{coefficienttype(src)}(os_grid)
    ArrayOperator(M, src, dest)
end

export levelweighed_approximate
function levelweighed_approximate(f, P::Platform, N::NTuple{D,Int}; options...) where D
    Nmin = minimum(N)
    @assert Nmin >= 2
    NN = ntuple(k->N[k]-Nmin+2,Val(D))
    if haskey(options,:L)
        q = div.(options[:L],1 .<< N)
        F,A,b,c,S = approximate(f, P, NN;options..., L=q.*(1 .<< NN))
    else
        F,A,b,c,S = approximate(f, P, NN;options...)
    end
    e = norm(A*c-b)
    w = [e]

    for i in NN[1]+1:N[1]
        NN = NN .+ 1
        W = LevelWeightingOperator(basis(dictionary(P,NN)),w)
        if haskey(options,:L)
            q = div.(options[:L],1 .<< N)
            F,A,b,c,S = approximate(f,P,NN;options..., L=q.*(1 .<< NN), weightedAZ=true, AZ_Cweight=W)
        else
            F,A,b,c,S = approximate(f,P,NN;options..., L=q.*(1 .<< NN), weightedAZ=true, AZ_Cweight=W)
        end
        e = norm(A*c-b)
        append!(w,e)
    end
    return F,A,b,c,S
end
end

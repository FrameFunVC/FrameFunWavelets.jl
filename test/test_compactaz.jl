module TestCompactAZ


using FrameFunWavelets, FrameFun, DomainSets, FrameFunTranslates, LinearAlgebra, Test, StaticArrays, SparseArrays
using WaveletsEvaluation.DWT: Wvl, Scl, Prl
using FrameFunWavelets.WaveletBases: GenericPeriodicEquispacedTranslates
@testset "nonzero_coefficients" begin
    P = CDFPlatform{3,3,Float64,Prl,Wvl,false}()
    N = 10
    dict = dictionary(P,N)
    g = sampling_grid(P,N;L=2048)
    A = AZ_A(P,N;L=2048);; Zt = AZ_Zt(P,N;L=2048);
    S = evaluation_operator(GenericPeriodicEquispacedTranslates(scalingbasis(dict)), g)
    iDWT = InverseDiscreteWaveletTransform(dict)
    @test A ≈ S*iDWT

    P = ExtensionFramePlatform(CDFPlatform{3,3,Float64,Prl,Wvl,true}(),0.0..0.4)
    N = 10
    g = sampling_grid(P,N;L=2048);
    dict = basis(dictionary(P,N))
    S = evaluation_operator(GenericPeriodicEquispacedTranslates(scalingbasis(dict)), g)
    iDWT = InverseDiscreteWaveletTransform(dict)
    A = AZ_A(P,N;L=2048); Z = AZ_Z(P,N;L=2048);
    @test A≈S*iDWT

    S2 = evaluation_operator(dest(basis(src(Z))),g)
    iDWT2 = inv(iDWT)'
    @test element(S2,1)'element(S,1)≈IdentityOperator(dict)
    @test Z≈S2*iDWT2
    @test element(Z,2)≈element(S2,1)

    spline_nonzero_coefs = nonzero_coefficients(ExtensionFramePlatform(BSplinePlatform{Float64,2}(),0.0..0.4),1<<N;L=2048)
    s = zeros(dict)
    s[spline_nonzero_coefs] .= NaN
    nonzero_coefs1 = findall(isnan.(iDWT'*s))
    spline_nonzero_coefs_refs = findall((sum(abs.(Matrix(S-S*S2'S));dims=1).>1e-12)[:])
    @test sort(spline_nonzero_coefs_refs)==sort(spline_nonzero_coefs)

    @test A-A*Z'A ≈ (S-S*S2'S)*iDWT
    M = Matrix(A-A*Z'*A)

    nonzero_coeffs = nonzero_coefficients(P,N;L=2048)
    MM = copy(M)
    MM[:,nonzero_coeffs].=0
    @test norm(MM) <1e-11

    P = ExtensionFramePlatform(NdCDFPlatform(2,3,3),(0.0..0.4)^2)
    N = (5,5);L=2 .*(1 .<<N)
    A = AZ_A(P,N;L=L); Zt = AZ_Zt(P,N;L=L);
    M = Matrix(A-A*Zt*A)
    nonzero_coeffs = nonzero_coefficients(P,N;L=L)
    MM = copy(M)
    MM[:,LinearIndices(1 .<< N)[nonzero_coeffs]].=0
    @test norm(MM) <1e-12

    P = ExtensionFramePlatform(NdCDFPlatform(2,4,4),(0.0..0.4)^2)
    N = (5,5);L=2 .*(1 .<<N)
    A = AZ_A(P,N;L=L); Zt = AZ_Zt(P,N;L=L);
    M = Matrix(A-A*Zt*A)
    nonzero_coeffs = nonzero_coefficients(P,N;L=L)
    MM = copy(M)
    MM[:,LinearIndices(1 .<< N)[nonzero_coeffs]].=0
    @test norm(MM) <1e-12

    P = ExtensionFramePlatform(DaubechiesPlatform(3),0.0..0.4)
    N = 4; L=2 .*(1 .<<N)
    A = AZ_A(P,N;L=L); Zt = AZ_Zt(P,N;L=L);
    M = Matrix(A-A*Zt*A)
    nonzero_coeffs = nonzero_coefficients(P,N;L=L)
    MM = copy(M)
    MM[:,nonzero_coeffs].=0
    @test norm(MM) <1e-10

    P = ExtensionFramePlatform(NdDaubechiesPlatform(2,3),(0.0..0.4)^2)
    N = (5,5);L=2 .*(1 .<<N)
    A = AZ_A(P,N;L=L); Zt = AZ_Zt(P,N;L=L);
    M = Matrix(A-A*Zt*A)
    nonzero_coeffs = nonzero_coefficients(P,N;L=L)
    MM = copy(M)
    MM[:,LinearIndices(1 .<< N)[nonzero_coeffs]].=0
    @test norm(MM) <1e-11
end

using FrameFunWavelets, FrameFun, Test, DomainSets, FrameFunTranslates, SparseArrays
@testset "reducedAAZAoperator" begin
    P = ExtensionFramePlatform(DaubechiesPlatform(2),0.0..0.4)
    N = 8
    M = Matrix(firstAZstepoperator(P,N;solverstyle=ReducedAZStyle()))
    rM = Matrix(reducedAAZAoperator(P,N;solverstyle=ReducedAZStyle()))
    sM = droptol!(sparse(rM),1e-12)
    srM = droptol!(sparse(M),1e-12)
    @test any(size(sM).< size(srM))
    @test nnz(sM)==nnz(srM)
    @test norm(M)≈norm(rM)

    P = ExtensionFramePlatform(NdCDFPlatform(2,3,3),(0.0..0.4)^2)
    N = 5,5
    M = Matrix(firstAZstepoperator(P,N;solverstyle=ReducedAZStyle()))
    rM = Matrix(reducedAAZAoperator(P,N;solverstyle=ReducedAZStyle()))
    sM = droptol!(sparse(rM),1e-12)
    srM = droptol!(sparse(M),1e-12)
    @test any(size(sM).< size(srM))
    @test nnz(sM)==nnz(srM)
    @test norm(M)≈norm(rM)


    P = ExtensionFramePlatform(NdCDFPlatform(2,4,4),(0.0..0.4)^2)
    N = 5,5
    M = Matrix(firstAZstepoperator(P,N;solverstyle=ReducedAZStyle()))
    rM = Matrix(reducedAAZAoperator(P,N;solverstyle=ReducedAZStyle()))
    sM = droptol!(sparse(rM),1e-12)
    srM = droptol!(sparse(M),1e-12)
    @test any(size(sM).< size(srM))
    @test nnz(sM)==nnz(srM)
    @test norm(M)≈norm(rM)
end

using FrameFunWavelets
using FrameFunWavelets.sparseiDWTEs: sparseiDWTMatrix, sparseiDWTE
@testset "2-D sparseiDWTE" begin
    dict = CDFWaveletBasis(4,4,2)^2
    sm = sparseiDWTMatrix(dict)
    M = Matrix(InverseDiscreteWaveletTransform(dict))
    @test M≈sm
    for i in eachindex(dict)
        @test sparseiDWTE(sm, [i])≈M[LinearIndices(size(dict))[i],:]'
    end
end

@testset "sparse_reducedAAZAoperator" begin
    P = ExtensionFramePlatform(DaubechiesPlatform{2,Float64,Wvl,false}(),0.0..0.4)
    N = 8
    M = Matrix(firstAZstepoperator(P,N))
    rM = sparse_reducedAAZAoperator(P,N;nz_tol=1e-12).A
    sM = droptol!(sparse(rM),1e-12)
    @test nnz(rM)==nnz(sM)
    srM = droptol!(sparse(M),1e-12)
    @test nnz(rM)==nnz(srM)
    @test (M)≈(rM)

    P = ExtensionFramePlatform(NdCDFPlatform(2,3,3),(0.0..0.4)^2)
    N = 5,5
    M = Matrix(firstAZstepoperator(P,N))
    rM = sparse_reducedAAZAoperator(P,N;nz_tol=1e-12).A
    sM = droptol!(sparse(rM),1e-12)
    @test nnz(rM)==nnz(sM)
    srM = droptol!(sparse(M),1e-12)
    @test nnz(rM)==nnz(srM)
    @test (M)≈(rM)


    P = ExtensionFramePlatform(CDFPlatform(4,4),(0.0..0.4))
    N = 5
    M = Matrix(firstAZstepoperator(P,N))
    rM = sparse_reducedAAZAoperator(P,N;nz_tol=1e-12).A
    sM = droptol!(sparse(rM),1e-12)
    @test nnz(rM)==nnz(sM)
    srM = droptol!(sparse(M),1e-12)
    @test norm(rM-Matrix(M),Inf)<1e-12
    @test nnz(rM)==nnz(srM)
    norm((M)-(rM))


    PBS = NdCDBSplinePlatform((3,3))
    PBSframe = ExtensionFramePlatform(PBS,(0.0..0.4)^2)
    Pbasis = NdCDFPlatform(2,4,4)
    P = ExtensionFramePlatform(Pbasis,(0.0..0.4)^2)
    N = 5,5
    L = 128,128
    dict = dictionary(P,N)
    ABS = AZ_A(PBSframe,1 .<<N; L=L)
    ZBS = AZ_Z(PBSframe,1 .<<N; L=L)
    A = AZ_A(P,N; L=L)
    Z = AZ_Z(P,N; L=L)

    iDWT = InverseDiscreteWaveletTransform(basis(dict))

    @test ABS*iDWT≈A
    @test ZBS*inv(iDWT)'≈Z
    @test A-A*Z'A≈firstAZstepoperator(P,N;L=L)

    M = Matrix(firstAZstepoperator(P,N;L=L))
    rM = sparse_reducedAAZAoperator(P,N;L=L,nz_tol=1e-12).A

    sM = droptol!(sparse(rM),1e-12)
    @test nnz(rM)==nnz(sM)
    srM = droptol!(sparse(M),1e-12)
    @test nnz(rM)==nnz(srM)
    @test norm((M)-(rM)) < 1e-11
end

D = .4*disk() + SVector(.5,.5)
    Pbasis = NdDaubechiesPlatform(2,2)
    P1 = ExtensionFramePlatform(Pbasis, D)
    N1 = 5,5
    m1 = (2,2)
    f1 = (x,y)->exp(x*y)

    D = .4*disk() + SVector(.5,.5)
    Pbasis = NdCDFPlatform(2,3,3)
    P3 = ExtensionFramePlatform(Pbasis, D)
    N3 = 5,5
    m3 = (2,2)
    f3 = (x,y)->exp(x*y)

    D = .4*disk() + SVector(.5,.5)
    Pbasis = NdCDFPlatform(2,4,2)
    P2 = ExtensionFramePlatform(Pbasis, D)
    N2 = 5,5
    m2 = (2,3)
    f2 = (x,y)->exp(x*y)

    D = .4*ball() + SVector(.5,.5,.5)
    Pbasis = NdDaubechiesPlatform(3,2)
    P4 = ExtensionFramePlatform(Pbasis, D)
    N4 = 4,4,3
    m4 = (2,2,2)
    f4 = (x,y,z)->exp(x*y*z)

    D = .4*ball() + SVector(.5,.5,.5)
    Pbasis = NdCDFPlatform(3,4,2)
    P5 = ExtensionFramePlatform(Pbasis, D)
    N5 = 3,3,3
    m5 = (2,2,2)
    f5 = (x,y,z)->exp(x*y*z)


@testset "sparse_reducedAAZAoperator" begin
    for (P,N,m) in zip((P1,P2,P3,P4,P5), (N1,N2,N3,N4,N5), (m1,m2,m3,m4,m5))
        op1 = firstAZstepoperator(P,N;L=m.*(1 .<< N))
        ops = sparse_reducedAAZAoperator(P,N;L=m.*(1 .<< N),nz_tol=1e-12)
        @test op1≈ops
        rM = ops.A
        sM = droptol!(sparse(rM),1e-12)
        @test nnz(rM)==nnz(sM)
    end
end


@testset "SparseAZ: FrameFunInterface" begin
    err = [1e-1,1e-4,1e-6,1e-1,2e-1]
    for (P,N,m,f,e) in zip((P1,P2,P3,P4,P5), (N1,N2,N3,N4,N5), (m1,m2,m3,m4,m5), (f1,f2,f3,f4,f5), err)
        F,A,b,c,_ = approximate(f, P, N; L=m.*(1 .<< N),solverstyle=SparseAZStyle())
        @test norm(A*c-b) < e
    end
end

end

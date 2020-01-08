
using FrameFunTranslates, FrameFunWavelets, DomainSets, StaticArrays, FrameFun, SparseArrays

module spmult
    using SparseArrays, LinearAlgebra
    using SparseArrays: prefer_sort

    function nnzspmatmul(A::SparseMatrixCSC, B::SparseMatrixCSC)
        mA, nA = size(A)
        nB = size(B, 2)
        nA == size(B, 1) || throw(DimensionMismatch())
        rowvalA = rowvals(A)
        rowvalB = rowvals(B)

        r = 0
        @inbounds begin
            ip = 1
            xb = fill(false, mA)
            nzindices = Vector{Int}(undef, mA)
            for i in 1:nB
                nzindex = 1
                k0 = ip - 1
                for jp in nzrange(B, i)
                    j = rowvalB[jp]
                    for kp in nzrange(A, j)
                        k = rowvalA[kp]
                        r = max(r, k+k0)
                        if !xb[k]
                            xb[k] = true
                            nzindices[nzindex] = k
                            nzindex += 1
                            ip += 1
                        end
                    end
                end
                for vp = 1:nzindex-1
                    xb[nzindices[vp]] = false
                end
            end

            return r
        end
    end

    function myspmatmul(A::SparseMatrixCSC{Tv,Ti}, B::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
        mA, nA = size(A)
        nB = size(B, 2)
        nA == size(B, 1) || throw(DimensionMismatch())

        rowvalA = rowvals(A); nzvalA = nonzeros(A)
        rowvalB = rowvals(B); nzvalB = nonzeros(B)
        @time nnzC = nnzspmatmul(A, B)

        colptrC = Vector{Ti}(undef, nB+1)
        rowvalC = Vector{Ti}(undef, nnzC)
        nzvalC = Vector{Tv}(undef, nnzC)
        nzpercol = nnzC ÷ max(nB, 1)

        @inbounds begin
            ip = 1
            xb = fill(false, mA)
            for i in 1:nB
                colptrC[i] = ip0 = ip
                k0 = ip - 1
                for jp in nzrange(B, i)
                    nzB = nzvalB[jp]
                    j = rowvalB[jp]
                    for kp in nzrange(A, j)
                        nzC = nzvalA[kp] * nzB
                        k = rowvalA[kp]
                        if xb[k]
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

        resize!(rowvalC, ip - 1)
        resize!(nzvalC, ip - 1)

        # This modification of Gustavson algorithm has sorted row indices
        C = SparseMatrixCSC(mA, nB, colptrC, rowvalC, nzvalC)
        return C
    end

end

using .spmult

L = 10
ns1 = round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1)))
ns2 = 3:10
ns3 = 2:6
    nsds = [ns1,ns2,ns3]

ds = 2:4

nnzA = zeros(Int,3, L,length(ds)).*NaN
    nnzAZ = copy(nnzA)
    nnzAZ1 = copy(nnzA)
    nnzAZtol = copy(nnzA)
    size1AZ1 = copy(nnzA)
    size2AZ1 = copy(nnzA)
    nnzAZ2 = copy(nnzA)
    size1AZ2 = copy(nnzA)
    size2AZ2 = copy(nnzA)
    size1A = copy(nnzA)
    size2A = copy(nnzA)
    include("fill_data_db_nnz.jl")
a

# You can leave this computation out
for (d,) in zip(1:3,), (i,n) in enumerate(nsds[d]), (j,p) in enumerate(ds[1:3]) #
    if d==1
        N = n
        Pbasis = DaubechiesPlatform(p)
    else
        N = ntuple(k->n,Val(d))
        Pbasis = NdDaubechiesPlatform(d, p)
    end
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    @show (1 .<< N), p

    if d==1
        Pbasis = DaubechiesPlatform(p)
        P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    elseif d==2
        Pbasis = NdDaubechiesPlatform(d, p)
        P = ExtensionFramePlatform(Pbasis, .4*disk() + SVector(.5,.5))
    elseif d==3
        Pbasis = NdDaubechiesPlatform(d, p)
        P = ExtensionFramePlatform(Pbasis, .4*ball() + SVector(.5,.5,.5))
    end

    A1, A2 = sparse_reducedAAZAoperator(P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver, return_scalingandidwteoperators=true)
    nnzAZ1[d,i,j] = nnz(A1)
    size1AZ1[d,i,j], size2AZ1[d,i,j] = size(A1)
    nnzAZ2[d,i,j] = nnz(A2')
    size1AZ2[d,i,j], size2AZ2[d,i,j] = size(A2)
    if (d==1 || (d==2 && (p<4 || (p==4 && n<10))) || d==3)
        M = spmult.myspmatmul(A1,copy(A2))
        nnzAZ[d,i,j] = nnz(M)
        nnzAZtol[d,i,j] = nnz(droptol!(M,1e-12))
    end
    if (d==1 || (d==2 && n < 10) || (d==3 && (n<5)))
       A = inv(solver(P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)).A
       nnzA[d,i,j] = nnz(A)
       size1A[d,i,j], size2A[d,i,j] = size(A)
    end

    @show nnzA[d,:,:]
    @show nnzAZ[d,:,:]
    @show nnzAZ1[d,:,:]
    @show nnzAZtol[d,:,:]
    @show size1AZ1[d,:,:]
    @show size2AZ1[d,:,:]
    @show nnzAZ2[d,:,:]
    @show size1AZ2[d,:,:]
    @show size2AZ2[d,:,:]
    @show size1A[d,:,:]
    @show size2A[d,:,:]
end

@show n1 = (1 .<< ns1)
@show n2 = (1 .<< ns2).^2
@show n3 = (1 .<< ns3).^3



for d in 1:3
    @show nnzA[d,:,:]
    @show nnzAZ[d,:,:]
    @show nnzAZ1[d,:,:]
    @show nnzAZtol[d,:,:]
    @show size1AZ1[d,:,:]
    @show size2AZ1[d,:,:]
    @show nnzAZ2[d,:,:]
    @show size1AZ2[d,:,:]
    @show size2AZ2[d,:,:]
    @show size1A[d,:,:]
    @show size2A[d,:,:]
end


using PGFPlotsX, LaTeXStrings, Printf, DocumentPGFPlots
ns = [n1,n2,n3]

PA = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],nnzA[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],nnzAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([ns[d],nnzA[d,:,p][1:length(ns[d])]]))
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],nnzAZtol[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([ns[d],nnzAZ[d,:,p][1:length(ns[d])]]))
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],nnzAZtol[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([ns[d],nnzA[d,:,p][1:length(ns[d])]]))
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)











fillA = nnzA./(size1A.*size2A)
fillAZ = nnzAZ./(size1AZ1.*size2AZ2)
fillAZtol = nnzAZtol./(size1AZ1.*size2AZ2)


PA = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],fillA[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],fillAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([ns[d],fillA[d,:,p][1:length(ns[d])]]))
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],fillAZtol[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([ns[d],fillAZ[d,:,p][1:length(ns[d])]]))
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],fillAZtol[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([ns[d],fillA[d,:,p][1:length(ns[d])]]))
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)



PA = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],nnzA[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
nnzAZ1
# Plot in function of number of splines that overlap with boundary
PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([size2AZ1[d,:,p][1:length(ns[d])],nnzAZ1[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    # vcat([
    # [Plot({style="black,dashed"},Table([size2AZ1[d,:,p][1:length(ns[d])],nnzA[d,:,p][1:length(ns[d])]]))
    # ] for p in 1:3
    # ]...)...,
    ] for d in 1:3]...)...)
# Plot in function of number of splines that overlap with boundary
PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([size2AZ1[d,:,p][1:length(ns[d])],nnzAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([size2AZ1[d,:,p][1:length(ns[d])],nnzAZ1[d,:,p][1:length(ns[d])]]))
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

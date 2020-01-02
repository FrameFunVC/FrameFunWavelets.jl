
using FrameFunTranslates, FrameFunWavelets, DomainSets, StaticArrays, FrameFun, SparseArrays

L = 10
ns1 = round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1)))
ns2 = 3:10
ns3 = 2:6
    nsds = [ns1,ns2,ns3]

ds = 2:4

nnzA = zeros(3, L,length(ds))
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
    M = A1*A2
    nnzAZ[d,i,j] = nnz(M)
    nnzAZtol[d,i,j] = nnz(droptol!(M,1e-12))

    A = inv(solver(P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)).A
    nnzA[d,i,j] = nnz(A)
    size1A[d,i,j], size2A[d,i,j] = size(A)
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

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],timingsAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)"))] for p in 1:3
    ]...)...,
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZ[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZ[d]))"))
    ] for d in 1:3]...)...)

PA = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZR1[d,1:length(ns[d]),p]])),
    LegendEntry("db"*latexstring("$(p+1)"))] for p in 1:3
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZR1[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),p]])) for p in 1:3]...,
    ] for d in 1:3]...)...)

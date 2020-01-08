
using FrameFunTranslates, FrameFunWavelets, DomainSets, StaticArrays, FrameFun, SparseArrays

L = 10
ns1 = round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1)))
ns2 = 3:10#round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2)))
ns3 = 2:6# round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3)))
    nsds = [ns1,ns2,ns3]

ds = ((2,4), (3,1), (3,5), (4,2))

f1 = (x)->exp(x)
f2 = (x,y)->exp(x*y)
f3 = (x,y,z)->exp(x*y*z)
fs = [f1,f2,f3]

timingsAZ = zeros(length(fs), L,length(ds))
    errorsAZ = copy(timingsAZ)
    allocAZ = copy(timingsAZ)
    timingsAZR1 = copy(timingsAZ)
    errorsAZR1 = copy(timingsAZ)
    allocAZR1 = copy(timingsAZ)
    timingsAZR2 = copy(timingsAZ)
    errorsAZR2 = copy(timingsAZ)
    allocAZR2 = copy(timingsAZ)
    timingsAZS = copy(timingsAZ)
    errorsAZS = copy(timingsAZ)
    allocAZS = copy(timingsAZ)
    timingsAZS2 = copy(timingsAZ)
    errorsAZS2 = copy(timingsAZ)
    allocAZS2 = copy(timingsAZ)
    timingsAZS3 = copy(timingsAZ)
    errorsAZS3 = copy(timingsAZ)
    allocAZS3 = copy(timingsAZ)

    include("fill_data_cdf.jl")

# You can leave this computation out
for (d,f) in zip(1:3,fs[1:3]), (i,n) in enumerate(nsds[d]), (j,t) in enumerate(ds[1:end]) #
    p,q = t
    if d==1
        N = n
        Pbasis = CDFPlatform(p,q)
    else
        N = ntuple(k->n,Val(d))
        Pbasis = NdCDFPlatform(d, p, q)
    end
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    @show (1 .<< N), p, q


    if (d==1) || (d==2 && n <= 8) || (d==3 && n<=5)

        F, timingsAZR1[d,i,j], allocAZR1[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false)
        errorsAZR1[d,i,j] = norm(F[2]*F[4]-F[3])
        if (i==1)
             F, timingsAZR1[d,i,j], allocAZR1[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=ReducedAZStyle())
        end

        @show timingsAZR1[d,:,:]
    end

    @show errorsAZR1
    @show allocAZR1

    if (d==1) || (d==2 && n <= 8) || (d==3 && n<=5)

        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=ReducedAZStyle(), directsolver=pSVD_solver, verbose=false)
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])
        if (i==1)
             F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=ReducedAZStyle(), directsolver=pQR_solver)
        end

        @show timingsAZR2[d,:,:]
    end

    @show errorsAZR2
    @show allocAZR2
    #
    if true
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])
        if (i==1)
             F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)
        end

        @show timingsAZS[d,:,:]
    end
    @show errorsAZS
    @show allocAZS

    if d==1
        Pbasis = CDFPlatform(p, q)
        P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    elseif d==2
        Pbasis = NdCDFPlatform(d, p, q)
        P = ExtensionFramePlatform(Pbasis, .4*disk() + SVector(.5,.5))
    elseif d==3
        Pbasis = NdCDFPlatform(d, p, q)
        P = ExtensionFramePlatform(Pbasis, .4*ball() + SVector(.5,.5,.5))
    end
    #
    if (d==1 || (d==2 && n<=9) || (d==3 && n <= 5))
        F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        errorsAZS2[d,i,j] = norm(F[2]*F[4]-F[3])
        if (i==1)
             F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)
        end

        @show timingsAZS2[d,:,:]
    end
    @show errorsAZS2
    @show allocAZS2

    if (d==1 || (d==2 && n<=7) || (d==3 && n <= 5))
        F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
        errorsAZS3[d,i,j] = norm(F[2]*F[4]-F[3])
        if (i==1)
             F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=DirectStyle(), directsolver=SPQR_solver)
        end

        @show timingsAZS3[d,:,:]
    end

    @show errorsAZS3
    @show allocAZS3

    if (d==1 && ((p,q) != (2,4) || ((p,q)==(2,4) && n <= 13))) ||
        (d==2 && (((p,q) == (2,4) && n<=6) || ((p,q) != (2,4) && n<=8))) ||
        (d==3 && n<=5)
        F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=AZStyle(), REG=pSVD_solver, verbose=false)
        errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])
        if (i==1)
            F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=2 .*(1 .<< N),threshold=1e-12,solverstyle=AZStyle(), REG=pQR_solver)
        end

        @show timingsAZ[d,:,:]
    end
    @show errorsAZ
    @show allocAZ

    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show timingsAZR1[d,:,:]
    @show errorsAZR1[d,:,:]
    @show allocAZR1[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show timingsAZS2[d,:,:]
    @show errorsAZS2[d,:,:]
    @show allocAZS2[d,:,:]
    @show timingsAZS3[d,:,:]
    @show errorsAZS3[d,:,:]
    @show allocAZS3[d,:,:]
end

@show n1 = (1 .<< ns1)
@show n2 = (1 .<< ns2).^2
@show n3 = (1 .<< ns3).^3
for d in 3:3
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show timingsAZR1[d,:,:]
    @show errorsAZR1[d,:,:]
    @show allocAZR1[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show timingsAZS[d,:,:]
    @show errorsAZS[d,:,:]
    @show allocAZS[d,:,:]
    @show timingsAZS2[d,:,:]
    @show errorsAZS2[d,:,:]
    @show allocAZS2[d,:,:]
    @show timingsAZS3[d,:,:]
    @show errorsAZS3[d,:,:]
    @show allocAZS3[d,:,:]
end

using PGFPlotsX, LaTeXStrings, Printf, DocumentPGFPlots
ns = [n1,n2,n3]

# Timings
asymptoticAZ = [
    (1e-6n1),
    (1e-6n2.^2),
    (1e-8n3.^((3*3-2)/3))
]
cAZ = map(x->@sprintf("N^{%1.2f}",x), [1,2,(3*3-2)/3]);#cAZ[1] *="\\log(N)"

PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],timingsAZ[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZ[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZ[d]))"))
    ] for d in 1:3]...)...)

asymptoticAZR1 = [
    (1e-6n1.*log.(ns1)),
    (1e-5n2.^(3/2)),
    (1e-7n3.^((3(3-1))/3))
]
cAZR1 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)]);cAZR1[1] *="\\log(N)"
PAZR1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZR1[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZR1[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)

asymptoticAZR2 = [
    1e-6n1,
    (5e-7n2.^(3/2))[1:8],
    1e-7n3.^((3(3-1))/3)
]
timingsAZR2[2,9,1] = 0
cAZR2 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)])
PAZR2  = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZR2[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZR1[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)

asymptoticAZS = [
    1e-5n1,
    1e-4n2,
    1e-4n3
]

cAZS = map(x->@sprintf("N^{%1.2f}",x), [1,1,1])
PAZS = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZS[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    # # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZS[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZS[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZR2[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)

# PAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],timingsAZS2[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", AZ")] for p in 1:3
#     ]...)...,
#     vcat([
#     [Plot(Table([ns[d],timingsAZS3[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)

imgpath = splitdir(@__FILE__())[1]








DocumentPGFPlots.savefigs(joinpath(imgpath,"AZtimings1d-3d_cdf"), PAZ)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR1timings1d-3d_cdf"), PAZR1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR2timings1d-3d_cdf"), PAZR2)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZStimings1d-3d_cdf"), PAZS)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSAStimings1d-3d_cdf"), PAZS1)



PAZI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZI")] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],timingsAZSI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
PeAZI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],errorsAZI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZI")] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],errorsAZSI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
PAZSI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", I")] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],timingsAZSI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
PeAZSI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],errorsI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", I")] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],errorsAZSI[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)


cc  = [1e-3,1e-3,1e-3]
@pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],errorsAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("p=$p")),
    ] for p in 1:3
    ]...)...,
    vcat([
    [Plot({style="black,dashed"},Table([ns[d],cc[p]*ns[d].^(Float64(-p))])),
    ] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
@pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],errorsAZR1[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
@pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],errorsAZR2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
@pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],errorsAZS[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
PeAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south west"},
    vcat([
    [Plot(Table([ns[d],errorsAZS2[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],errorsAZS3[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)
# DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSASerrors1d-3d"), PeAZS1)


using FrameFunTranslates, FrameFunWavelets, DomainSets, StaticArrays, FrameFun, SparseArrays

L = 10
ns1 = round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1)))
ns2 = 3:10#round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2)))
ns3 = 2:6# round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3)))
    nsds = [ns1,ns2,ns3]

ds = 2:4

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

    mAZR2 = zeros(Int,size(timingsAZ)...)
    nAZR2 = copy(mAZR2)
    rAZR2 = copy(mAZR2)
    MAZR2 = copy(mAZR2)
    NAZR2 = copy(mAZR2)
    r1 = copy(timingsAZ)
    r2 = copy(timingsAZ)
    include("fill_data_db.jl")

# You can leave this computation out
for (d,f) in zip(1:3,fs[1:3]), (i,n) in enumerate(nsds[d]), (j,p) in zip(3:3,ds[3:3]) #
    if d==1
        N = n
        Pbasis = DaubechiesPlatform(p)
        L = 4 * (1 << N)
    else
        N = ntuple(k->n,Val(d))
        Pbasis = NdDaubechiesPlatform(d, p)
        L = 4 .* (1 .<< N)
    end
    threshold = p == 4 ? 1e-4 : 1e-12
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    @show (1 .<< N), p, threshold
    S = true_nonzero_reducedAAZAoperator(P,N;L=L)
    @show r1[d,i,j] = size(S.operators[end], 1)/size(S.operators[1], 2)

    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=AZStyle(), directsolver=pQR_solver, verbose=true)
        end
        F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=AZStyle(), directsolver=pQR_solver, verbose=false)
        errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZ[d,:,:]
        @show errorsAZ[d,:,:]
        @show allocAZ[d,:,:]
    end

    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), directsolver=pQR_solver, verbose=true)
        end
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), directsolver=pQR_solver, verbose=false)
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZR2[d,:,:]
        @show errorsAZR2[d,:,:]
        @show allocAZR2[d,:,:]
    end


    if (d==1 || (d==2 && n<=9) || (d==3 && n <= 5))
        if (i==1)
             F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=SparseAZStyle(), directsolver=SPQR_solver,verbose=true)
        end
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=SparseAZStyle(), verbose=false, directsolver=SPQR_solver)
        errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])


        @show timingsAZS[d,:,:]
        @show errorsAZS[d,:,:]
        @show allocAZS[d,:,:]
    end

    if d==1
        Pbasis = DaubechiesPlatform(p)
        P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
        L = 4 * (1<<N)
    elseif d==2
        Pbasis = NdDaubechiesPlatform(d, p)
        P = ExtensionFramePlatform(Pbasis, .4*disk() + SVector(.5,.5))
        L = 2 .* (1 .<<N)
    elseif d==3
        Pbasis = NdDaubechiesPlatform(d, p)
        P = ExtensionFramePlatform(Pbasis, .4*ball() + SVector(.5,.5,.5))
        L = 2 .* (1 .<<N)
    end
    S = true_nonzero_reducedAAZAoperator(P,N;L=L)
    @show r2[d,i,j] = size(S.operators[end], 1)/size(S.operators[1], 2)

    if (d==1 || (d==2 && n<=9) || (d==3 && n <= 5))
        if (i==1)
             F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=SparseAZStyle(), directsolver=SPQR_solver,verbose=true)
        end
        F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=SparseAZStyle(), verbose=false, directsolver=SPQR_solver)
        errorsAZS2[d,i,j] = norm(F[2]*F[4]-F[3])


        @show timingsAZS2[d,:,:]
        @show errorsAZS2[d,:,:]
        @show allocAZS2[d,:,:]
    end

    if ((d==1 && n<=16) || (d==2 && n<=6) || (d==3 && n <= 5))
        if (i==1)
             F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=DirectStyle(), directsolver=SPQR_solver, verbose=true)
        end
        F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
        errorsAZS3[d,i,j] = norm(F[2]*F[4]-F[3])


        @show timingsAZS3[d,:,:]
        @show errorsAZS3[d,:,:]
        @show allocAZS3[d,:,:]
    end

    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
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
    println()
end

for (d,f) in zip(1:3,fs[1:3]), (i,n) in enumerate(nsds[d]), (j,p) in enumerate(ds[1:end]) #
    if d==1
        N = n
        Pbasis = DaubechiesPlatform(p)
        L = 4 * (1 << N)
    else
        N = ntuple(k->n,Val(d))
        Pbasis = NdDaubechiesPlatform(d, p)
        L = 4 .* (1 .<< N)
    end
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    @show (1 .<< N), p
    if (d==1) || (d==2 && n <= 8) || (d==3 && n<=5)
        S = solver(P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), directsolver=pSVD_solver)
        a = S.psolver.operators[3]
        @show mAZR2[d,i,j], nAZR2[d,i,j], rAZR2[d,i,j] = size(a)..., size(a.solver.S,1)
        @show MAZR2[d,i,j], NAZR2[d,i,j]= size(AZ_A(P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), directsolver=pSVD_solver))
    end
end

@show n1 = (1 .<< ns1)
@show n2 = (1 .<< ns2).^2
@show n3 = (1 .<< ns3).^3
for d in 1:3
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
    @show rAZR2[d,:,:]
    @show mAZR2[d,:,:]
    @show nAZR2[d,:,:]
    @show MAZR2[d,:,:]
    @show NAZR2[d,:,:]
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
    [Plot(Table([ns[d],timingsAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry("db"*latexstring("$(p+1)"))] for p in 1:3
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
    [Plot(Table([ns[d],timingsAZR1[d,1:length(ns[d]),p]])),
    LegendEntry("db"*latexstring("$(p+1)"))] for p in 1:3
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZR1[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),p]])) for p in 1:3]...,
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
    [Plot(Table([ns[d],timingsAZR2[d,1:length(ns[d]),p]])),
    LegendEntry("db"*latexstring("$(p+1)"))] for p in 1:3
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZR2[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),p]])) for p in 1:3]...,
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
    [Plot(Table([ns[d],timingsAZS[d,1:length(ns[d]),p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    # # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d],asymptoticAZS[d][1:length(ns[d])]])),
    LegendEntry(latexstring("\\mathcal O($(cAZS[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZR2[d,1:length(ns[d]),p]])) for p in 1:3]...,
    ] for d in 1:3]...)...)

PAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZS2[d,:,p]])),
    LegendEntry(latexstring("p=$p")*", AZ")] for p in 1:3
    ]...)...,
    vcat([
    [Plot(Table([ns[d],timingsAZS3[d,:,p]])),
    LegendEntry(latexstring("p=$p"))] for p in 1:3
    ]...)...,
    ] for d in 1:3]...)...)

imgpath = splitdir(@__FILE__())[1]
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZtimings1d-3d"), PAZ)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR1timings1d-3d"), PAZR1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR2timings1d-3d"), PAZR2)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZStimings1d-3d"), PAZS)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSAStimings1d-3d"), PAZS1)



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

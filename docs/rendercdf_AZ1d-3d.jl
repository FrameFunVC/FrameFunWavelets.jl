
using FrameFunTranslates, FrameFunWavelets, DomainSets, StaticArrays, FrameFun, SparseArrays

L = 10
ns1 = round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1)))
ns2 = 3:10#round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2)))
ns3 = 2:6# round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3)))
    nsds = [ns1,ns2,ns3]

ds = ((3,1), (3,5), (4,2))

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

    include("fill_data_cdf.jl")

# You can leave this computation out
for (d,f) in zip(1:3,fs[1:3]), (i,n) in enumerate(nsds[d]), (j,t) in enumerate(ds[1:end]) #
    p,q = t
    if d==1
        N = n
        Pbasis = CDFPlatform(p,q)
        L = 4 * (1 << N)
    else
        N = ntuple(k->n,Val(d))
        Pbasis = NdCDFPlatform(d, p, q)
        L = 4 .* (1 .<< N)
    end
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    @show (1 .<< N), p, q
    S = true_nonzero_reducedAAZAoperator(P,N;L=L)
    @show r1[d,i,j] = size(S.operators[end], 1)/size(S.operators[1], 2)

    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=AZStyle(), directsolver=pQR_solver, verbose=true)
        end
        F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=AZStyle(), directsolver=pQR_solver, verbose=false)
        errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZ[d,:,:]
        @show errorsAZ[d,:,:]
        @show allocAZ[d,:,:]
    end

    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=ReducedAZStyle(), directsolver=pQR_solver, verbose=true)
        end
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=ReducedAZStyle(), directsolver=pQR_solver, verbose=false)
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])

        @show timingsAZR2[d,:,:]
        @show errorsAZR2[d,:,:]
        @show allocAZR2[d,:,:]
    end


    if (d==1 || (d==2 && n<=9) || (d==3 && n <= 5))
        if (i==1)
             F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=SparseAZStyle(), directsolver=SPQR_solver,verbose=true)
        end
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, directsolver=SPQR_solver)
        errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])


        @show timingsAZS[d,:,:]
        @show errorsAZS[d,:,:]
        @show allocAZS[d,:,:]
    end

    if d==1
        Pbasis = CDFPlatform(p, q)
        P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
        L = 4 * (1<<N)
    elseif d==2
        Pbasis = NdCDFPlatform(d, p, q)
        P = ExtensionFramePlatform(Pbasis, .4*disk() + SVector(.5,.5))
        L = 2 .* (1 .<<N)
    elseif d==3
        Pbasis = NdCDFPlatform(d, p, q)
        P = ExtensionFramePlatform(Pbasis, .4*ball() + SVector(.5,.5,.5))
        L = 2 .* (1 .<<N)
    end
    S = true_nonzero_reducedAAZAoperator(P,N;L=L)
    @show r2[d,i,j] = size(S.operators[end], 1)/size(S.operators[1], 2)

    if (d==1 || (d==2 && n<=9) || (d==3 && n <= 5))
        if (i==1)
             F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=SparseAZStyle(), directsolver=SPQR_solver,verbose=true)
        end
        F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, directsolver=SPQR_solver)
        errorsAZS2[d,i,j] = norm(F[2]*F[4]-F[3])


        @show timingsAZS2[d,:,:]
        @show errorsAZS2[d,:,:]
        @show allocAZS2[d,:,:]
    end

    if ((d==1 && n<=16) || (d==2 && n<=6) || (d==3 && n <= 5))
        if (i==1)
             F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=DirectStyle(), directsolver=SPQR_solver, verbose=true)
        end
        F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
        errorsAZS3[d,i,j] = norm(F[2]*F[4]-F[3])


        @show timingsAZS3[d,:,:]
        @show errorsAZS3[d,:,:]
        @show allocAZS3[d,:,:]
    end

    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
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



# for (d,f) in zip(1:3,fs[1:3]), (i,n) in enumerate(nsds[d]), (j,t) in enumerate(ds[1:end]) #
#     p,q = t
#     if d==1
#         N = n
#         Pbasis = CDFPlatform(p,q)
#     else
#         N = ntuple(k->n,Val(d))
#         Pbasis = NdCDFPlatform(d, p, q)
#     end
#     P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
#     @show (1 .<< N), p, q
#     if (d==1) || (d==2 && n <= 8) || (d==3 && n<=5)
#         S = solver(P, N;L=L,threshold=1e-12,solverstyle=ReducedAZStyle(), directsolver=pSVD_solver)
#         a = S.psolver.operators[3]
#         @show mAZR2[d,i,j], nAZR2[d,i,j], rAZR2[d,i,j] = size(a)..., size(a.solver.S,1)
#         @show MAZR2[d,i,j], NAZR2[d,i,j]= size(AZ_A(P, N;L=L,threshold=1e-12,solverstyle=ReducedAZStyle(), directsolver=pSVD_solver))
#     end
# end

@show n1 = (1 .<< ns1)
@show n2 = (1 .<< ns2).^2
@show n3 = (1 .<< ns3).^3
for d in 1:3
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
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
    (1e-5n1),
    (1e-7n2[1:6].^2),
    (1e-7n3[1:4].^((3*3-2)/3))
]
cAZ = map(x->@sprintf("N^{%1.2f}",x), [1,4/2,(3*3-2)/3]);
function myplot(d)
    if d==3
        return @pgf [Plot({style="black,dotted,thick"},Table([ns[d][1:length(asymptoticAZ[d])],3asymptoticAZ[d].^(2/2.33)])),
                    LegendEntry(latexstring("\\mathcal O($(@sprintf("N^{%1.2f}",2)))"))]
    elseif d==2
        return @pgf [Plot({style="black,dotted,thick"},Table([ns[d][1:length(asymptoticAZ[d])],3asymptoticAZ[d].^(1.5/2)])),
                    LegendEntry(latexstring("\\mathcal O($(@sprintf("N^{%1.2f}",1.5)))"))]
    else
        return []
    end
end
PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    vcat([
    [Plot(Table([ns[d],timingsAZ[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZ[d])],asymptoticAZ[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZ[d]))")),
    myplot(d)...
    ] for d in 1:3]...,
    )...)

asymptoticAZR1 = [
    (1e-6n1.*log.(ns1)),
    (1e-5n2.^(3/2).*log.(ns2).^2)[1:6],
    (1e-6n3.^((3(3-1))/3).*log.(ns3).^2)[1:4]
]
cAZR1 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)]);cAZR1[1] *="\\log(N)";cAZR1[2] *="\\log(N)^2";cAZR1[3] *="\\log(N)^2"
PAZR1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZR1[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZR1[d])],asymptoticAZR1[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)

asymptoticAZR2 = [
    1e-5n1,
    (5e-6n2.^(3/2))[1:6],
    (1e-7n3.^((3(3-1))/3))[1:4]
]

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
    Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZR2[d])],asymptoticAZR2[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR2[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)

asymptoticAZS = [
    1e-5n1,
    1e-4n2[1:7],
    1e-3n3[1:4]
]
timingsAZS[3,5,1] = 0
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
        Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZS[d])],asymptoticAZS[d]])),
        LegendEntry(latexstring("\\mathcal O($(cAZS[d]))")),
        # Plot previous plot
        [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
        ] for d in 1:3]...)...)


PAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZS3[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))
    ] for (k,pq) in enumerate(ds)
    ]...)...,
    vcat([
    [Plot({style="black,dotted"},Table([ns[d],timingsAZS2[d,:,k][1:length(ns[d])]])),
    # LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2]), AZ"))
    ]    for (k,pq) in enumerate(ds)
    ]...)...,
    ] for d in 1:3]...)...)



imgpath = splitdir(@__FILE__())[1]

DocumentPGFPlots.savefigs(joinpath(imgpath,"AZtimings1d-3d_cdf"), PAZ)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR1timings1d-3d_cdf"), PAZR1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR2timings1d-3d_cdf"), PAZR2)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZStimings1d-3d_cdf"), PAZS)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSAStimings1d-3d_cdf"), PAZS1)


s = 1
Pn = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
    vcat([
    [PlotInc(Table([ns[d],rAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry("r")] for (k,pq) in zip(s:s,ds[s:s])
    ]...)...,
    vcat([
    [PlotInc(Table([ns[d],nAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry("n")] for (k,pq) in zip(s:s,ds[s:s])
    ]...)...,
    vcat([
    [PlotInc(Table([ns[d],mAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry("m")] for (k,pq) in zip(s:s,ds[s:s])
    ]...)...,
    vcat([
    [PlotInc(Table([ns[d],NAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry("N")] for (k,pq) in zip(s:s,ds[s:s])
    ]...)...,
    vcat([
    [PlotInc(Table([ns[d],MAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry("M")] for (k,pq) in zip(s:s,ds[s:s])
    ]...)...,
    ] for d in 1:3]...)...)


# Errors
PAZ = @pgf GroupPlot({scale_only_axis, width=".5\\textwidth",height=".35\\textwidth",group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",ymin=1e-15,ymax=1e0,legend_cell_align="left",legend_pos="north east"},
    vcat([
    [PlotInc({},Table([ns[d],errorsAZ[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    ] for d in 1:3]...,
    )...)
PAZR1 = @pgf GroupPlot({scale_only_axis, width=".5\\textwidth",height=".35\\textwidth",group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",ymin=1e-15,ymax=1e0,legend_cell_align="left",legend_pos="north east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],errorsAZR1[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)
PAZR2  = @pgf GroupPlot({scale_only_axis, width=".5\\textwidth",height=".35\\textwidth",group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",ymin=1e-15,ymax=1e0,legend_cell_align="left",legend_pos="north east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],errorsAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)
PAZS = @pgf GroupPlot({scale_only_axis, width=".5\\textwidth",height=".35\\textwidth",group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",ymin=1e-15,ymax=1e0,legend_cell_align="left",legend_pos="north east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],errorsAZS[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
    ] for d in 1:3]...)...)
PAZS1 = @pgf GroupPlot({scale_only_axis, width=".5\\textwidth",height=".35\\textwidth",group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],errorsAZS3[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))
    ] for (k,pq) in enumerate(ds)
    ]...)...,
    vcat([
    [Plot({style="black,dotted"},Table([ns[d],errorsAZS2[d,:,k][1:length(ns[d])]])),
    # LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2]), AZ"))
    ]    for (k,pq) in enumerate(ds)
    ]...)...,
    ] for d in 1:3]...)...)

imgpath = splitdir(@__FILE__())[1]

DocumentPGFPlots.savefigs(joinpath(imgpath,"AZerrors1d-3d_cdf"), PAZ)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR1errors1d-3d_cdf"), PAZR1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR2errors1d-3d_cdf"), PAZR2)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSerrors1d-3d_cdf"), PAZS)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSASerrors1d-3d_cdf"), PAZS1)



PAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
    vcat([[
    {xlabel=latexstring("N"),xmode="log",ymode="log",ymin=1e-15,ymax=1e0,legend_cell_align="left",legend_pos="north east"},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],errorsAZS2[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2]), AZ"))] for (k,pq) in enumerate(ds)
    ]...)...,
    vcat([
    [Plot(Table([ns[d],errorsAZS3[d,:,k][1:length(ns[d])]])),
    LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
    ]...)...,
    ] for d in 1:3]...)...)


# PAZI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],timingsAZI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", AZI")] for p in 1:3
#     ]...)...,
#     vcat([
#     [Plot(Table([ns[d],timingsAZSI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# PeAZI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],errorsAZI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", AZI")] for p in 1:3
#     ]...)...,
#     vcat([
#     [Plot(Table([ns[d],errorsAZSI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# PAZSI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],timingsI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", I")] for p in 1:3
#     ]...)...,
#     vcat([
#     [Plot(Table([ns[d],timingsAZSI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# PeAZSI = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="north west"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],errorsI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", I")] for p in 1:3
#     ]...)...,
#     vcat([
#     [Plot(Table([ns[d],errorsAZSI[d,:,p]])),
#     LegendEntry(latexstring("p=$p")*", AZSI")] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)


# cc  = [1e-3,1e-3,1e-3]
# @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     vcat([
#     [Plot(Table([ns[d],errorsAZ[d,:,p][1:length(ns[d])]])),
#     LegendEntry(latexstring("p=$p")),
#     ] for p in 1:3
#     ]...)...,
#     vcat([
#     [Plot({style="black,dashed"},Table([ns[d],cc[p]*ns[d].^(Float64(-p))])),
#     ] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     vcat([
#     [Plot(Table([ns[d],errorsAZR1[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     vcat([
#     [Plot(Table([ns[d],errorsAZR2[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     vcat([
#     [Plot(Table([ns[d],errorsAZS[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# PeAZS1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south west"},
#     vcat([
#     [Plot(Table([ns[d],errorsAZS2[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     vcat([
#     [Plot(Table([ns[d],errorsAZS3[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     ] for d in 1:3]...)...)
# DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSASerrors1d-3d"), PeAZS1)

# asymptoticAZR1 = [
#     (1e-6n1.*log.(ns1)),
#     (1e-5n2.^(3/2).*log.(ns2).^2)[1:6],
#     (1e-6n3.^((3(3-1))/3).*log.(ns3).^2)[1:4]
# ]
# cAZR1 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)]);cAZR1[1] *="\\log(N)";cAZR1[2] *="\\log(N)^2";cAZR1[3] *="\\log(N)^2"
# PAZR1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],timingsAZR1[d,:,k][1:length(ns[d])]])),
#     LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
#     ]...)...,
#     # Add asymptotic lines
#     Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZR1[d])],asymptoticAZR1[d]])),
#     LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
#     # Plot previous plot
#     [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
#     ] for d in 1:3]...)...)

# asymptoticAZR2 = [
#     1e-6n1,
#     (5e-7n2.^(3/2))[1:8],
#     1e-7n3.^((3(3-1))/3)
# ]
#
# cAZR2 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)])
# PAZR2  = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],timingsAZR2[d,:,k][1:length(ns[d])]])),
#     LegendEntry("cdf"*latexstring("$(pq[1])$(pq[2])"))] for (k,pq) in enumerate(ds)
#     ]...)...,
#     # Add asymptotic lines
#     Plot({style="black,dashed"},Table([ns[d],asymptoticAZR2[d][1:length(ns[d])]])),
#     LegendEntry(latexstring("\\mathcal O($(cAZR2[d]))")),
#     # Plot previous plot
#     [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ds)]...,
#     ] for d in 1:3]...)...)

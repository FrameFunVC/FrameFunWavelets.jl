using FrameFunTranslates, FrameFunWavelets, DomainSets, StaticArrays, FrameFun, SparseArrays, WaveletsEvaluation, LowRankApprox

L = 10
ns1 = round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1)))
ns2 = 3:10#round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2)))
ns3 = 2:5# round.(Int,log2.((10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3)))
    nsds = [ns1,ns2,ns3]

ws = (db2,db3,db4,cdf31,cdf33,cdf35,cdf42,cdf44,cdf46)
Ds = ((0.0..0.5),.35disk()+SVector(.5,.5),.4ball()+SVector(.5,.5,.5))


f1 = (x)->exp(x)
f2 = (x,y)->exp(x*y)
f3 = (x,y,z)->exp(x*y*z)
fs = [f1,f2,f3]

timingsAZ = zeros(length(fs), L,length(ws))
    errorsAZ = copy(timingsAZ)
    allocAZ = copy(timingsAZ)
    coefnormAZ = copy(timingsAZ)
    timingsAZsvd = copy(timingsAZ)
    errorsAZsvd = copy(timingsAZ)
    allocAZsvd = copy(timingsAZ)
    coefnormAZsvd = copy(timingsAZ)
    timingsAZR1 = copy(timingsAZ)
    errorsAZR1 = copy(timingsAZ)
    allocAZR1 = copy(timingsAZ)
    coefnormAZR1 = copy(timingsAZ)
    timingsAZR2 = copy(timingsAZ)
    errorsAZR2 = copy(timingsAZ)
    allocAZR2 = copy(timingsAZ)
    coefnormAZR2 = copy(timingsAZ)
    timingsAZR3 = copy(timingsAZ)
    errorsAZR3 = copy(timingsAZ)
    allocAZR3 = copy(timingsAZ)
    coefnormAZR3 = copy(timingsAZ)
    timingsAZS = copy(timingsAZ)
    errorsAZS = copy(timingsAZ)
    allocAZS = copy(timingsAZ)
    coefnormAZS = copy(timingsAZ)
    timingsS = copy(timingsAZ)
    errorsS = copy(timingsAZ)
    allocS = copy(timingsAZ)
    coefnormS = copy(timingsAZ)

    threshold=1e-12
    lraoptions=LRAOptions(atol=1e-14,rtol=1e-12)
    include("fill_data.jl")

# You can leave this computation out
for (d,f) in zip(1:3,fs[1:3]), (i,n) in enumerate(nsds[d]), (j,w) in enumerate(ws[1:end]) #
    if d==1
        N = n
        Pbasis = waveletplatform(w)
        L = 2 * (1 << N)
    else
        N = ntuple(k->n,Val(d))
        Pbasis = ProductPlatform(ntuple(k->waveletplatform(w), Val(d))...)
        L = 2 .* (1 .<< N)
    end
    P = ExtensionFramePlatform(Pbasis, Ds[d])
    @show (1 .<< N), w

    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=AZStyle(), REG=pQR_solver, lraoptions=lraoptions, verbose=false)
        end
        F, timingsAZ[d,i,j], allocAZ[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=AZStyle(), REG=pQR_solver, lraoptions=lraoptions, verbose=false)
        errorsAZ[d,i,j] = norm(F[2]*F[4]-F[3])
        coefnormAZ[d,i,j] = norm(F[4])

        @show timingsAZ[d,:,:]
        @show errorsAZ[d,:,:]
    end

    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZsvd[d,i,j], allocAZsvd[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=AZStyle(), REG=pSVD_solver, lraoptions=lraoptions, verbose=false)
        end
        F, timingsAZsvd[d,i,j], allocAZsvd[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=AZStyle(), REG=pSVD_solver, lraoptions=lraoptions, verbose=false)
        errorsAZsvd[d,i,j] = norm(F[2]*F[4]-F[3])
        coefnormAZsvd[d,i,j] = norm(F[4])

        @show timingsAZsvd[d,:,:]
        @show errorsAZsvd[d,:,:]
    end


    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZR1[d,i,j], allocAZR1[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), REG=QR_solver, verbose=false)
        end
        F, timingsAZR1[d,i,j], allocAZR1[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), REG=QR_solver, verbose=false)
        errorsAZR1[d,i,j] = norm(F[2]*F[4]-F[3])
        coefnormAZR1[d,i,j] = norm(F[4])

        @show timingsAZR1[d,:,:]
        @show errorsAZR1[d,:,:]
    end


    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), REG=pQR_solver, lraoptions=lraoptions, verbose=false)
        end
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), REG=pQR_solver, lraoptions=lraoptions, verbose=false)
        errorsAZR2[d,i,j] = norm(F[2]*F[4]-F[3])
        coefnormAZR2[d,i,j] = norm(F[4])

        @show timingsAZR2[d,:,:]
        @show errorsAZR2[d,:,:]
    end

    if (d==1) ||
        (d==2 && n<=8) ||
        (d==3 && n<=5)
        if (i==1)
            F, timingsAZR3[d,i,j], allocAZR3[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), REG=pSVD_solver, lraoptions=lraoptions, verbose=false)
        end
        F, timingsAZR3[d,i,j], allocAZR3[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), REG=pSVD_solver, lraoptions=lraoptions, verbose=false)
        errorsAZR3[d,i,j] = norm(F[2]*F[4]-F[3])
        coefnormAZR3[d,i,j] = norm(F[4])

        @show timingsAZR3[d,:,:]
        @show errorsAZR3[d,:,:]
    end


    if (d==1 || (d==2 && n<=9) || (d==3 && n <= 5))
        if (i==1)
             F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=SparseAZStyle(), REG=SPQR_solver,verbose=false)
        end
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        errorsAZS[d,i,j] = norm(F[2]*F[4]-F[3])
        coefnormAZS[d,i,j] = norm(F[4])


        @show timingsAZS[d,:,:]
        @show errorsAZS[d,:,:]
    end

    if ((d==1 && n<=16) || (d==2 && n<=6) || (d==3 && n <= 5))
        if (i==1)
             F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=DirectStyle(), directsolver=SPQR_solver, verbose=false)
        end
        F, timingsS[d,i,j], allocS[d,i,j], _ = @timed approximate(f, P, N;L=L,threshold=threshold,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
        errorsS[d,i,j] = norm(F[2]*F[4]-F[3])
        coefnormS[d,i,j] = norm(F[4])


        @show timingsS[d,:,:]
        @show errorsS[d,:,:]
    end

    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show coefnormAZ[d,:,:]
    @show timingsAZsvd[d,:,:]
    @show errorsAZsvd[d,:,:]
    @show allocAZsvd[d,:,:]
    @show coefnormAZsvd[d,:,:]
    # @show timingsAZR1[d,:,:]
    # @show errorsAZR1[d,:,:]
    # @show allocAZR1[d,:,:]
    # @show coefnormAZR1[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show coefnormAZR2[d,:,:]
    @show timingsAZR3[d,:,:]
    @show errorsAZR3[d,:,:]
    @show allocAZR3[d,:,:]
    @show coefnormAZR3[d,:,:]
    # @show timingsAZS[d,:,:]
    # @show errorsAZS[d,:,:]
    # @show allocAZS[d,:,:]
    # @show coefnormAZS[d,:,:]
    # @show timingsS[d,:,:]
    # @show errorsS[d,:,:]
    # @show allocS[d,:,:]
    # @show coefnormS[d,:,:]
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
#         S = solver(P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), directsolver=pSVD_solver)
#         a = S.psolver.operators[3]
#         @show mAZR2[d,i,j], nAZR2[d,i,j], rAZR2[d,i,j] = size(a)..., size(a.solver.S,1)
#         @show MAZR2[d,i,j], NAZR2[d,i,j]= size(AZ_A(P, N;L=L,threshold=threshold,solverstyle=ReducedAZStyle(), directsolver=pSVD_solver))
#     end
# end

@show n1 = (1 .<< ns1)
@show n2 = (1 .<< ns2).^2
@show n3 = (1 .<< ns3).^3
for d in 1:3
    println()
    @show timingsAZ[d,:,:]
    @show errorsAZ[d,:,:]
    @show allocAZ[d,:,:]
    @show coefnormAZ[d,:,:]
    @show timingsAZsvd[d,:,:]
    @show errorsAZsvd[d,:,:]
    @show allocAZsvd[d,:,:]
    @show coefnormAZsvd[d,:,:]
    # @show timingsAZR1[d,:,:]
    # @show errorsAZR1[d,:,:]
    # @show allocAZR1[d,:,:]
    # @show coefnormAZR1[d,:,:]
    @show timingsAZR2[d,:,:]
    @show errorsAZR2[d,:,:]
    @show allocAZR2[d,:,:]
    @show coefnormAZR2[d,:,:]
    @show timingsAZR3[d,:,:]
    @show errorsAZR3[d,:,:]
    @show allocAZR3[d,:,:]
    @show coefnormAZR3[d,:,:]
    # @show timingsAZS[d,:,:]
    # @show errorsAZS[d,:,:]
    # @show allocAZS[d,:,:]
    # @show coefnormAZS[d,:,:]
    # @show timingsS[d,:,:]
    # @show errorsS[d,:,:]
    # @show allocS[d,:,:]
    # @show coefnormS[d,:,:]
    println()
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
PAZ = @pgf GroupPlot({
    legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="north west",
    xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",
    group_style={group_size={"3 by 1"}}},
    vcat([[
    {
        # (d==3 ? {} : {legend_to_name="test"*string(randn()),})...,
        },
    vcat([
    [PlotInc(Table([ns[d],timingsAZ[d,:,k][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (k,pq) in enumerate(ws)
    ]...)...,
    Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZ[d])],asymptoticAZ[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZ[d]))")),
    myplot(d)...,
    ] for d in 1:3]...,
    )...)

asymptoticAZR1 = [
    (1e-6n1.*log.(ns1)),
    (1e-5n2.^(3/2).*log.(ns2).^2)[1:6],
    (1e-6n3.^((3(3-1))/3).*log.(ns3).^2)[1:4]
]
cAZR1 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)]);cAZR1[1] *="\\log(N)";cAZR1[2] *="\\log(N)^2";cAZR1[3] *="\\log(N)^2"
PAZR1 = @pgf GroupPlot({
    legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="north west",
    xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",
    group_style={group_size={"3 by 1"}},
    },
    vcat([[
    {
    # (d==3 ? {} : {legend_to_name="test"*string(randn()),})...,
    },
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZR1[d,:,k][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (k,pq) in enumerate(ws)
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZR1[d])],asymptoticAZR1[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ws)]...,
    ] for d in 1:3]...)...)

asymptoticAZR2 = [
    1e-5n1,
    (5e-6n2.^(3/2))[1:6],
    (1e-7n3.^((3(3-1))/3))[1:4]
]

cAZR2 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)])
PAZR2  = @pgf GroupPlot({
        legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="north west",
        xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",
        group_style={group_size={"3 by 1"}},
    },
    vcat([[
    {},
    # Plot this plot
    vcat([
    [Plot(Table([ns[d],timingsAZR2[d,:,k][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (k,pq) in enumerate(ws)
    ]...)...,
    # Add asymptotic lines
    Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZR2[d])],asymptoticAZR2[d]])),
    LegendEntry(latexstring("\\mathcal O($(cAZR2[d]))")),
    # Plot previous plot
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ws)]...,
    ] for d in 1:3]...)...)
asymptoticAZS = [
    1e-5n1,
    1e-4n2[1:7],
    1e-3n3[1:4]
]
cAZS = map(x->@sprintf("N^{%1.2f}",x), [1,1,1])
PAZS = @pgf GroupPlot({
        legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="north west",
        xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",
        group_style={group_size={"3 by 1"}},
    },
        vcat([[
        {},
        # Plot this plot
        vcat([
        [Plot(Table([ns[d],timingsAZS[d,:,k][1:length(ns[d])]])),
        LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (k,pq) in enumerate(ws)
        ]...)...,
        # # Add asymptotic lines
        Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZS[d])],asymptoticAZS[d]])),
        LegendEntry(latexstring("\\mathcal O($(cAZS[d]))")),
        # Plot previous plot
        [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ws)]...,
        ] for d in 1:3]...)...)
PS = @pgf GroupPlot({
        legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="north west",
        xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",
        group_style={group_size={"3 by 1"}},
    },
        vcat([[
        {},
        # Plot this plot
        vcat([
        [Plot(Table([ns[d],timingsS[d,:,k][1:length(ns[d])]])),
        LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (k,pq) in enumerate(ws)
        ]...)...,
        # # Add asymptotic lines
        Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZS[d])],asymptoticAZS[d]])),
        LegendEntry(latexstring("\\mathcal O($(cAZS[d]))")),
        # Plot previous plot
        [Plot({style="black,dotted"},Table([ns[d],timingsAZS[d,1:length(ns[d]),k]])) for (k,pq) in enumerate(ws)]...,
        ] for d in 1:3]...)...)

Ptime = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
    legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="north west",ymin=1e-4,ymax=1e4,
    width="0.5\\textwidth",height="0.4\\textwidth",
    xlabel=L"N",ylabel="Time (s.)",
    group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",
        group_size={"3 by 5"}}},
    vcat([[
    {(d!=1 ? {legend_to_name="test"*string(randn())} : {})...,},
    vcat([
    [Plot(Table([ns[d],timingsAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],timingsAZR1[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],timingsAZR2[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],timingsAZS[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],timingsS[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],timingsAZS[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    )
# errorsAZS db4 does not converge
Presidual = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
    legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="south west",ymin=1e-15,ymax=10,
    width="0.5\\textwidth",height="0.4\\textwidth",
    xlabel=L"N",ylabel="Residual",
    group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",
        group_size={"3 by 5"}}},
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],errorsAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],errorsAZR1[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],errorsAZR2[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],errorsAZS[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {(d!=3 ? {legend_to_name="test"*string(randn())} : {})...,},
    vcat([
    [Plot(Table([ns[d],errorsS[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],errorsAZS[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    )

[support_length(Primal,wavelet,w) for w in ws]



Pcoefnorm = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
    legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="south west",
    width="0.5\\textwidth",height="0.4\\textwidth",
    # xlabel=L"N",ylabel=L"\|x\|",
    group_style={x_descriptions_at="edge bottom",
        # y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",
        group_size={"3 by 5"}}},
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],coefnormAZ[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],coefnormAZR1[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],coefnormAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],coefnormAZR2[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],coefnormAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([ns[d],coefnormAZS[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],coefnormAZ[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {(d!=3 ? {legend_to_name="test"*string(randn())} : {})...,},
    vcat([
    [Plot(Table([ns[d],coefnormS[d,:,p][1:length(ns[d])]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([ns[d],coefnormAZS[d,:,p][1:length(ns[d])]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    )




Pefficiency = @pgf GroupPlot({xmode="log",ymode="log",legend_cell_align="left",
    legend_style={fill_opacity=0.6,text_opacity=1},legend_pos="south west",ymin=1e-15,ymax=10,
    width="0.5\\textwidth",height="0.4\\textwidth",
    xlabel="Time (s.)",ylabel="Residual",
    group_style={x_descriptions_at="edge bottom",
        y_descriptions_at="edge left", horizontal_sep="1em",vertical_sep="1em",
        group_size={"3 by 5"}}},
    vcat([[
    {(d!=3 ? {legend_to_name="test"*string(randn())} : {})...,},
    vcat([
    [Plot(Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([timingsAZR1[d,:,p],errorsAZR1[d,:,p]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([timingsAZR2[d,:,p],errorsAZR2[d,:,p]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([timingsAZS[d,:,p],errorsAZS[d,:,p]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZ[d,:,p],errorsAZ[d,:,p]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    vcat([[
    {legend_to_name="test"*string(randn())},
    vcat([
    [Plot(Table([timingsS[d,:,p],errorsS[d,:,p]])),
    LegendEntry(latexstring("\\texttt{$(WaveletsEvaluation.DWT.name(pq))}"))] for (p,pq) in enumerate(ws)
    ]...)...,
    [Plot({style="black,dotted"},Table([timingsAZS[d,:,p],errorsAZS[d,:,p]])) for (p,pq) in enumerate(ws)]...
    ] for d in 1:3]...)...,
    )

imgpath = splitdir(@__FILE__())[1]

DocumentPGFPlots.savefigs(joinpath(imgpath,"AZtimings1d-3d"), PAZ)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR1timings1d-3d"), PAZR1)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR2timings1d-3d"), PAZR2)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZStimings1d-3d"), PAZS)
DocumentPGFPlots.savefigs(joinpath(imgpath,"AZAStimings1d-3d"), PS)

DocumentPGFPlots.savefigs(joinpath(imgpath,"timings1d-3d"), Ptime)
DocumentPGFPlots.savefigs(joinpath(imgpath,"errors1d-3d"), Presidual)
DocumentPGFPlots.savefigs(joinpath(imgpath,"efficiency1d-3d"), Pefficiency)
DocumentPGFPlots.savefigs(joinpath(imgpath,"coefnorms1d-3d"), Pcoefnorm)

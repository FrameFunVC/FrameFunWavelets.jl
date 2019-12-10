
using FrameFunTranslates, FrameFunWavelets, DomainSets, StaticArrays, FrameFun

L = 10
ns1 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1))
ns2 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e6))...,L)).^(1/2))
ns3 = round.(Int,(10.0.^LinRange(log10.((1e2, 1e5))...,L)).^(1/3))
    nsds = [ns1,ns2,ns3]

ds = 1:3

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

    timingsAZI = copy(timingsAZ)
    errorsAZI = copy(timingsAZ)
    allocAZI = copy(timingsAZ)
    timingsAZSI = copy(timingsAZ)
    errorsAZSI = copy(timingsAZ)
    allocAZSI = copy(timingsAZ)
    timingsI = copy(timingsAZ)
    errorsI = copy(timingsAZ)
    allocI = copy(timingsAZ)
    # include("fill_data.jl")

# You can leave this computation out
for (d,f) in zip(1:3,fs[1:3]), (i,n) in zip(1:10,nsds[d][1:10]), (j,p) in enumerate(ds[1:3]) #
    if d==1
        Pbasis = DaubechiesPlatform(p)
    else
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
    end
    P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)

    if d ==1 ||
           (d==2 &&((p==2 || p==3) && n <= 28) || ((p==1) &&n^d <= 5e5)) ||
           (d==3 &&((p==1 && n^d <=2e5) || (p==2 && n<=36) || (p==3 && n<=36)))

        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle(), verbose=false)
        errorsAZR2[d,i,j] = residual(f, F)
        F, timingsAZR2[d,i,j], allocAZR2[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=ReducedAZStyle())

        @show timingsAZR2
    end
    # @show errorsAZR2
    # @show allocAZR2

    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36))) #32
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        errorsAZS[d,i,j] = residual(f, F)
        F, timingsAZS[d,i,j], allocAZS[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)

        @show timingsAZS
    end
    # @show errorsAZS
    # @show allocAZS

    if d==1
        Pbasis = CDBSplinePlatform(p)
        P = ExtensionFramePlatform(Pbasis, (0.0..0.5)^d)
    elseif d==2
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
        P = ExtensionFramePlatform(Pbasis, .4*disk() + SVector(.5,.5))
    elseif d==3
        Pbasis = NdCDBSplinePlatform(ntuple(k->p,Val(d)))
        P = ExtensionFramePlatform(Pbasis, .4*ball() + SVector(.5,.5,.5))
    end

    if (d==1) || (d==2) || (d==3 && (p==1||p==2||(p==3&&n<=36)))
        F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), verbose=false, REG=SPQR_solver)
        errorsAZS2[d,i,j] = residual(f, F)
        F, timingsAZS2[d,i,j], allocAZS2[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=SparseAZStyle(), REG=SPQR_solver)

        @show timingsAZS2
    end
    # @show errorsAZS2
    # @show allocAZS2


    F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), verbose=false, directsolver=SPQR_solver)
    errorsAZS3[d,i,j] = residual(f, F)
    F, timingsAZS3[d,i,j], allocAZS3[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-12,solverstyle=DirectStyle(), directsolver=SPQR_solver)

    @show timingsAZS3

    # @show errorsAZS3
    # @show allocAZS3


    # F, timingsAZI[d,i,j], allocAZI[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-3,solverstyle=AZStyle(), verbose=false, REG=LSQR_solver)
    # errorsAZI[d,i,j] = residual(f, F)
    # F, timingsAZI[d,i,j], allocAZI[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-3,solverstyle=AZStyle(), REG=LSQR_solver)
    #
    # @show timingsAZI

    # @show errorsAZI
    # @show allocAZI

    F, timingsAZSI[d,i,j], allocAZSI[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-3,solverstyle=SparseAZStyle(), verbose=false, REG=LSQR_solver)
    errorsAZSI[d,i,j] = residual(f, F)
    F, timingsAZSI[d,i,j], allocAZSI[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-3,solverstyle=SparseAZStyle(), REG=LSQR_solver)

    @show timingsAZSI

    # @show errorsAZSI
    # @show allocAZSI

    F, timingsI[d,i,j], allocI[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-3,solverstyle=DirectStyle(), verbose=false, directsolver=LSQR_solver)
    errorsI[d,i,j] = residual(f, F)
    F, timingsI[d,i,j], allocI[d,i,j], _ = @timed Fun(f, P, N;threshold=1e-3,solverstyle=DirectStyle(), directsolver=LSQR_solver)

    @show timingsI

    # @show errorsI
    # @show allocI
end

@show n1 = ns1
@show n2 = ns2.^2
@show n3 = ns3.^3
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
    @show timingsAZI[d,:,:]
    @show errorsAZI[d,:,:]
    @show allocAZI[d,:,:]
    @show timingsI[d,:,:]
    @show errorsI[d,:,:]
    @show allocI[d,:,:]
    @show timingsAZSI[d,:,:]
    @show errorsAZSI[d,:,:]
    @show allocAZSI[d,:,:]
end

# using PGFPlotsX, LaTeXStrings, Printf, DocumentPGFPlots
# ns = [n1,n2,n3]
#
# # Timings
# asymptoticAZ = [
#     (1e-6n1.*log.(n1))[3:end],
#     (1e-6n2.^2)[1:6],
#     (1e-8n3.^((3*3-2)/3))[1:8]
# ]
# cAZ = map(x->@sprintf("N^{%1.2f}",x), [1,2,(3*3-2)/3]);cAZ[1] *="\\log(N)"
#
# PAZ = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     vcat([
#     [Plot(Table([ns[d][d==1 ? (3:end) : (:)],timingsAZ[d,d==1 ? (3:end) : (:),p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     Plot({style="black,dashed"},Table([ns[d][d==1 ? (3:end) : (:)][1:length(asymptoticAZ[d])],asymptoticAZ[d]])),
#     LegendEntry(latexstring("\\mathcal O($(cAZ[d]))"))
#     ] for d in 1:3]...)...)
#
# asymptoticAZR1 = [
#     (1e-6n1.*log.(ns1))[3:end],
#     (1e-5n2.^(3/2))[1:7],
#     (1e-7n3.^((3(3-1))/3))[1:9]
# ]
# cAZR1 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)]);cAZR1[1] *="\\log(N)"
# timingsAZR1[3,10,1] = 0
# PAZR1 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d][d==1 ? (3:end) : (:)],timingsAZR1[d,d==1 ? (3:end) : (:),p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     # Add asymptotic lines
#     Plot({style="black,dashed"},Table([ns[d][d==1 ? (3:end) : (:)][1:length(asymptoticAZR1[d])],asymptoticAZR1[d]])),
#     LegendEntry(latexstring("\\mathcal O($(cAZR1[d]))")),
#     # Plot previous plot
#     [Plot({style="black,dotted"},Table([ns[d][d==1 ? (3:end) : (:)],timingsAZ[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...,
#     ] for d in 1:3]...)...)
#
# asymptoticAZR2 = [
#     1e-6n1,
#     (5e-7n2.^(3/2))[1:8],
#     1e-7n3.^((3(3-1))/3)
# ]
# timingsAZR2[2,9,1] = 0
# cAZR2 = map(x->@sprintf("N^{%1.2f}",x), [1,3/2,((3(3-1))/3)])
# PAZR2 = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],timingsAZR2[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     # Add asymptotic lines
#     Plot({style="black,dashed"},Table([ns[d][1:length(asymptoticAZR2[d])],asymptoticAZR2[d]])),
#     LegendEntry(latexstring("\\mathcal O($(cAZR2[d]))")),
#     # Plot previous plot
#     [Plot({style="black,dotted"},Table([ns[d][d==1 ? (3:end) : (:)],timingsAZR1[d,d==1 ? (3:end) : (:),p]])) for p in 1:3]...,
#     ] for d in 1:3]...)...)
#
# asymptoticAZS = [
#     1e-5n1,
#     1e-4n2,
#     1e-4n3
# ]
#
# cAZS = map(x->@sprintf("N^{%1.2f}",x), [1,1,1])
# PAZS = @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xlabel=latexstring("N"),xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     # Plot this plot
#     vcat([
#     [Plot(Table([ns[d],timingsAZS[d,:,p]])),
#     LegendEntry(latexstring("p=$p"))] for p in 1:3
#     ]...)...,
#     # # Add asymptotic lines
#     Plot({style="black,dashed"},Table([ns[d],asymptoticAZS[d]])),
#     LegendEntry(latexstring("\\mathcal O($(cAZS[d]))")),
#     # Plot previous plot
#     [Plot({style="black,dotted"},Table([ns[d],timingsAZR2[d,:,p]])) for p in 1:3]...,
#     ] for d in 1:3]...)...)
#
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
#
# # imgpath = splitdir(@__FILE__())[1]
# # DocumentPGFPlots.savefigs(joinpath(imgpath,"AZtimings1d-3d"), PAZ)
# # DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR1timings1d-3d"), PAZR1)
# # DocumentPGFPlots.savefigs(joinpath(imgpath,"AZR2timings1d-3d"), PAZR2)
# # DocumentPGFPlots.savefigs(joinpath(imgpath,"AZStimings1d-3d"), PAZS)
# # DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSAStimings1d-3d"), PAZS1)
#
#
#
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
#
#
# cc  = [1e-3,1e-3,1e-3]
# @pgf GroupPlot({group_style={group_size={"3 by 1"}}},
#     vcat([[
#     {xmode="log",ymode="log",legend_cell_align="left",legend_pos="south east"},
#     vcat([
#     [Plot(Table([ns[d],errorsAZ[d,:,p]])),
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
# # DocumentPGFPlots.savefigs(joinpath(imgpath,"AZSASerrors1d-3d"), PeAZS1)

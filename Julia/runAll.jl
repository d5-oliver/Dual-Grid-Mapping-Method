using Plots
using LaTeXStrings

include("FVM.jl")
include("IM.jl")
include("DM.jl")
include("par.jl")

function runAll()

    x0::Float64 = 0
    xL::Float64 = 30

    N::Int64 = 101
    M::Int64 = 11

    T::Float64 = 18
    K::Int64 = 3600

    omega::Int64 = 1
    sigma::Float64 = 1e+10

    testCase::Int64 = 6

    plotFigs::Bool = true

    parVals::Any = constructPar(x0,xL,T,N,M,K,omega,sigma,testCase)

    c_FVM::Array{Float64} = FVM(parVals)
    c_IM::Array{Float64}, C_IM::Array{Float64} = IM(parVals)
    c_DM::Array{Float64}, C_DM::Array{Float64} = DM(parVals)

    if plotFigs

        gr(size=(800,600), xtickfontsize=13, ytickfontsize=13, xguidefontsize=20, yguidefontsize=24, legendfontsize=16, dpi=100)

        plot(parVals.xF, c_FVM[:,parVals.tSteps], lw=3, lc=:gray, label="", legend=false)
        plot!(parVals.xF,c_IM[:,parVals.tSteps], lw=3, lc=:red, ls=:dash, label="")
        xlabel!(L"x")
        ylabel!(L"$c_i^{(k)}, \widetilde{c}_i^{(k)}, C_m^{(k)}$")
        display(scatter!(parVals.xC, C_IM[:,parVals.tSteps], lw=3, markerstrokewidth=3, markerstrokecolor=:red, mc=:red, lc=:red, ls=:dot, label=""))

        plot(parVals.xF,c_FVM[:,parVals.tSteps], lw=3, lc=:gray, label="", legend=false)
        plot!(parVals.xF,c_DM[:,parVals.tSteps], lw=3, lc=:blue, ls=:dash, label="")
        xlabel!(L"x")
        ylabel!(L"$c_i^{(k)}, \widetilde{c}_i^{(k)}, C_m^{(k)}$")
        display(scatter!(parVals.xC, C_DM[:,parVals.tSteps], lw=3, markerstrokewidth=3, markerstrokecolor=:blue, mc=:blue, lc=:blue, ls=:dot, label=""))

        local plotInterval = 0:0.05:0.125
        local numEntries = length(plotInterval)

        plot(plotInterval, repeat([1],numEntries), lw=3, lc=:gray, label="", showaxis=false, grid=false, legend=false, xlims=(0,1), ylims=(-0.25,1.25))
        annotate!(0.3, 1, "Fine-Grid solution")

        plot!(plotInterval, repeat([0.75],numEntries), lw=3, lc=:red, ls=:dash, label="", grid=false, legend=false)
        annotate!(0.3, 0.75, "Fine-Grid IM solution")

        scatter!([0.05], [0.5], mode="markers", name="", markerstyle=:dot, markerstrokewidth=3, markercolor=:red, markerstrokecolor=:red)
        annotate!(0.3, 0.5, "Coarse-Grid IM solution")

        plot!(plotInterval, repeat([0.25],numEntries), lw=3, lc=:blue, label="", grid=false, legend=false)
        annotate!(0.3, 0.25, "Fine-Grid DM solution")

        scatter!([0.05], [0], mode="markers", name="", markerstyle=:dot, markerstrokewidth=3, markercolor=:blue, markerstrokecolor=:blue)
        display(annotate!(0.3, 0, "Coarse-Grid DM solution"))

    end

end

runAll()

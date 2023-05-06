using Plots
using LaTeXStrings

include("FVM.jl")
include("IM.jl")
include("DM.jl")
include("layers.jl")
include("par.jl")

# Execute the runAll() function to generate results for the fine-grid solution method and the dual-grid mapping method (IM and DM).
# Edit values here and in par.jl to modify test problems.
# Set testCase = 1, 2, ..., 6 to pick a test problem.
# Set plotFigs = true to plot the results.
# Please note: Plot output is saved to the current directory.

function runAll()

	testCase = 5

	layers = setLayers(testCase)

	x0 = 0.0
	xL = 30.0

	N = 3601
	xF = x0:(xL-x0)/(N-1):xL
	# If test case is a layered problem, ensure a node is placed on interfaces
	if testCase in [1,2,3,4]
		while !isempty(setdiff(layers[2:end-1,1],xF))
			N += 1
			xF = x0:(xL-x0)/(N-1):xL
		end
	end

	M = 11
	# Ensure coarse grid is a subset of fine grid nodes.
	while mod( (N-1)/(M-1) , 1 ) > 0
		M += 1;
	end
	xC = x0:(xL-x0)/(M-1):xL
	
	T = 18.0
	K = 3600

	omega = 1.0
	sigma = 1.0e+10

	plotFigs = true

	parVals = constructPar(x0,xL,T,N,xF,M,xC,K,omega,sigma,testCase,layers)

	c_FVM = FVM(parVals)
	c_IM, C_IM = IM(parVals)
	c_DM, C_DM = DM(parVals)

	if plotFigs

		gr(size=(800,600), xtickfontsize=13, ytickfontsize=13, xguidefontsize=20, yguidefontsize=24, legendfontsize=16, dpi=100)

		plot(parVals.xF, c_FVM[:,parVals.tSteps], lw=3, lc=:gray, label="", legend=false)
		plot!(parVals.xF,c_IM[:,parVals.tSteps], lw=3, lc=:red, ls=:dash, label="")
		xlabel!(L"x")
		ylabel!(L"$c_i^{(k)}, \widetilde{c}_i^{(k)}, C_m^{(k)}$")
		scatter!(parVals.xC, C_IM[:,parVals.tSteps], lw=3, markerstrokewidth=3, markerstrokecolor=:red, mc=:red, lc=:red, ls=:dot, label="")

		savefig("IM_Comparison.png")

		plot(parVals.xF,c_FVM[:,parVals.tSteps], lw=3, lc=:gray, label="", legend=false)
		plot!(parVals.xF,c_DM[:,parVals.tSteps], lw=3, lc=:blue, ls=:dash, label="")
		xlabel!(L"x")
		ylabel!(L"$c_i^{(k)}, \widetilde{c}_i^{(k)}, C_m^{(k)}$")
		scatter!(parVals.xC, C_DM[:,parVals.tSteps], lw=3, markerstrokewidth=3, markerstrokecolor=:blue, mc=:blue, lc=:blue, ls=:dot, label="")

		savefig("DM_Comparison.png")

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
		annotate!(0.3, 0, "Coarse-Grid DM solution")

		savefig("Legend.png")

	end

end

runAll()

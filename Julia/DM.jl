using LinearAlgebra

include("fvmStd.jl")
include("fvmSig.jl")

function constructG(M,P,totalNodes,tau,itMatLocal,dm)

	# Construct the mapping matrix G [Equation 25]
	Gm::Array{Float64} = zeros(totalNodes,M)

	idxStore::StepRangeLen{Int64} = 0:0

	for m in 1:M-1

		idxStore = (m-1)*P+1:m*P

		Gm[idxStore,[m,m+1]] = itMatLocal[idxStore,idxStore] \ (tau * [dm[1:P,m] dm[2:P+1,m]])

	end

	return Gm

end

function constructPhi(M,P,tau,itMatLocal,bm,c,phim)

	# Construct the phi vector [Equation 26]
	idxGrid::StepRangeLen{Int64} = 0:0
	idxStore::StepRangeLen{Int64} = 0:0

	for m in 1:M-1

		idxGrid = (m-1)*(P-1)+1:m*(P-1)+1
		idxStore = (m-1)*P+1:m*P

		phim[idxStore] = itMatLocal[idxStore,idxStore] \ (c[idxGrid] .+ tau*bm[:,m])

	end

	return phim

end

function DM(par)

	# Transport coefficients and problem parameters
	x0::Float64 = par.x0
	xL::Float64 = par.xL

	N::Int64 = par.N
	M::Int64 = par.M
	P::Int64 = (N-1)/(M-1) + 1
	totalNodes::Int64 = N + M - 2

	xF::StepRangeLen{Float64} = par.xF
	xC::StepRangeLen{Float64} = par.xC
	
	tau::Float64 = par.tau
	t::StepRangeLen{Float64} = par.t
	K::Int64 = par.K

	R::Any = par.R
	D::Any = par.D
	v::Any = par.v
	mu::Any = par.mu
	Gamma::Any = par.Gamma

	omega::Int64 = par.omega

	sigma::Float64 = par.sigma

	alpha0::Float64 = par.alpha0
	beta0::Float64 = par.beta0
	sigma0::Float64 = par.sigma0
	g0::Any = par.g0(t)

	alphaL::Float64 = par.alphaL
	betaL::Float64 = par.betaL
	sigmaL::Float64 = par.sigmaL
	gL::Any = par.gL(t)

	# Construct the iteration matrix A [Equation 7] and vector b [Equation 8]
	A::Array{Float64}, b::Array{Float64} = fvmStd(xF,R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL)
	itMatFull::Array{Float64} = I - tau*A

	# Construct the iteration matrices A [Equation 21] and vectors b [Equation 23] corresponding to locally defined sub-problems [Equations 16-19]
	Am::Array{Float64}, bm::Array{Float64}, dm::Array{Float64} = constructArrays(R,D,v,mu,Gamma,omega,sigma,xF,xC)
	itMatLocal::Array{Float64} = I - tau*Am

	b0::Float64 = b[1]
	bL::Float64 = b[N]

	b[2:N-1] = tau*b[2:N-1]

	val0::Float64 = Gamma(x0)/R(x0)
	valL::Float64 = Gamma(xL)/R(xL)

	idx::Array{Int64} = setdiff(1:totalNodes, (1:M-2)*P.+1)

	# Construct the mapping matrix G [Equation 25]
	Gm::Array{Float64} = constructG(M,P,totalNodes,tau,itMatLocal,dm)
	Gm = Gm[idx,:]
	GmT::Array{Float64} = transpose(Gm)

	itMatCoarse::Array{Float64} = GmT * itMatFull * Gm

	c::Array{Float64} = zeros(N,K+1)
	c[:,1] = par.c0(xF)

	C::Array{Float64} = zeros(M,K+1)
	C[:,1] = par.c0(xC)

	phim::Array{Float64} = zeros(totalNodes)

	# Time stepping
	for k in 1:K

		# Construct the phi vector [Equation 26]
		phim = constructPhi(M,P,tau,itMatLocal,bm,c[:,k],phim)

		b[1] = tau * (b0*g0[k+1] + val0)
		b[N] = tau * (bL*gL[k+1] + valL)

		# Approximate solution on the coarse grid [Equation 12]
		C[:,k+1] = itMatCoarse \ (GmT * (c[:,k] - itMatFull*phim[idx] + b))

		# Reconstruct solution on the fine grid [Equation 11]
		c[:,k+1] = Gm * C[:,k+1] + phim[idx]

	end

	return c, C

end
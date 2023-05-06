using LinearAlgebra

include("fvmStd.jl")
include("fvmSig.jl")

# Create and store the LU factors of the locally defined iteration matrices (which are defined in [Equation 20])
function luLocal(M,P,blockMat)

	itMatLU = Array{LU{Float64}}(undef,M-1,1)

	@views @inbounds for m in 1:M-1

		idxStore = (m-1)*P+1:m*P

		itMatLU[m] = lu(blockMat[idxStore,idxStore])

	end

	return itMatLU

end

function DM(par)

	# Transport coefficients and problem parameters
	x0 = par.x0
	xL = par.xL

	N = par.N
	M = par.M
	P = div(N-1,M-1) + 1
	totalNodes = N + M - 2

	xF = par.xF
	xC = par.xC

	tau = par.tau
	t = par.t
	K = par.K

	R = par.R
	D = par.D
	v = par.v
	mu = par.mu
	Gamma = par.Gamma

	omega = par.omega

	sigma = par.sigma

	alpha0 = par.alpha0
	beta0 = par.beta0
	sigma0 = par.sigma0
	g0 = par.g0(t)

	alphaL = par.alphaL
	betaL = par.betaL
	sigmaL = par.sigmaL
	gL = par.gL(t)

	# Construct the iteration matrix A [Equation 7] and vector b [Equation 8]
	A, b = fvmStd(xF,R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL)
	itMatFull = I - tau*A

	# Pre-allocate the b vectors for use in the time-stepping loop
	bval = zeros(N,K+1)
	bval[1,:] .= tau .* (b[1] .* g0 .+ Gamma(x0) / R(x0))
	bval[N,:] .= tau .* (b[N] .* gL .+ Gamma(xL) / R(xL))
	@views @inbounds for k in 1:K
		bval[2:N-1,k] .= tau .* b[2:N-1]
	end

	# Construct the iteration matrices A [Equation 21] and vectors b, u, and v [Equation 23] corresponding to locally defined sub-problems [Equations 16-19]
	Am, bm, uvm = constructArrays(R,D,v,mu,Gamma,omega,sigma,xF,xC)
	itMatLocal = I - tau*Am
	itMatLocal = luLocal(M,P,itMatLocal)

	# Pre-allocate bm and uvm arrays for use with constructing mapping matrix G and vectors phi
	bmVec = tau * bm
	uvmVec = tau * uvm

	# Set indices for the double-ups
	idx = setdiff(1:totalNodes, (1:M-2)*P.+1)

	# Pre-allocate temp arrays for use in constructing G and phi
	rhsTmpGm = zeros(P,2)
	rhsTmpPhim = zeros(P)

	# Construct the mapping matrix G [Equation 25]
	Gm = zeros(totalNodes,M)
	@views @inbounds for m in 1:M-1

		idxStore = (m-1)*P+1:m*P

		rhsTmpGm .= uvmVec[:,:,m]
		ldiv!(itMatLocal[m],rhsTmpGm)

		Gm[idxStore,[m,m+1]] .= rhsTmpGm

	end
	Gm = Gm[idx,:]
	GmT = transpose(Gm)

	# Construct the coarse-grid iteration matrix
	itMatCoarse = GmT * itMatFull * Gm
	itMatCoarse = lu(itMatCoarse)

	# Pre-allocate full Nx1 solutions
	c = zeros(N,K+1)
	c[:,1] = par.c0(xF)

	# Pre-allocate coarse Mx1 solutions
	C = zeros(M,K+1)
	C[:,1] = par.c0(xC)

	# Pre-allocate phi vector
	phim = zeros(totalNodes)

	# Pre-allocate RHS vectors for the full (stepTmpbFull) NxN problem, and the coarse (stepTmpbCoarse) MxM problem
	stepTmpbFull = similar(xF)
	stepTmpbCoarse = similar(xC)

	# Time stepping
	@views @inbounds for k in 1:K

		# Construct the phi vector [Equation 26]
		for m in 1:M-1

			idxGrid = (m-1)*(P-1)+1:m*(P-1)+1
			idxStore = (m-1)*P+1:m*P

			rhsTmpPhim .= c[idxGrid,k]
			rhsTmpPhim .+= bmVec[:,m]
			ldiv!(itMatLocal[m],rhsTmpPhim)

			phim[idxStore] .= rhsTmpPhim

		end

		# Construct the full Nx1 RHS vector
		stepTmpbFull .= c[:,k] 
		stepTmpbFull .-= itMatFull * phim[idx]
		stepTmpbFull .+= bval[:,k+1]

		# Construct the coarse Mx1 RHS vector
		stepTmpbCoarse .= GmT * stepTmpbFull

		# Approximate solution on the coarse grid [Equation 12]
		ldiv!(itMatCoarse,stepTmpbCoarse)
		C[:,k+1] .= stepTmpbCoarse

		# Reconstruct solution on the fine grid [Equation 11]
		stepTmpbFull .= Gm * C[:,k+1]
		stepTmpbFull .+= phim[idx]

		c[:,k+1] .= stepTmpbFull

	end

	return c, C

end
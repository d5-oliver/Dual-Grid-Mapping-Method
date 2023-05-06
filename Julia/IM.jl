using LinearAlgebra

include("fvmStd.jl")

# Construct mapping matrix by linear interpolation [Equation 13]
function interpG(N,xF,M,xC)

	G = zeros(N,M)
	P = div(N-1,M-1) + 1

	for m in 1:M-1

		idx = (m-1)*(P-1)+1:m*(P-1)+1

		G[idx,m] .= (xC[m+1] .- xF[idx]) ./ (xC[m+1] .- xC[m])
		G[idx,m+1] .= (xF[idx] .- xC[m]) ./ (xC[m+1] .- xC[m])

	end

	return G

end

function IM(par)

	# Transport coefficients and problem parameters
	x0 = par.x0
	xL = par.xL

	N = par.N
	M = par.M

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

	alpha0 = par.alpha0
	beta0 = par.beta0
	sigma0 = par.sigma0
	g0 = par.g0(t)

	alphaL = par.alphaL
	betaL = par.betaL
	sigmaL = par.sigmaL
	gL = par.gL(t)

	# Construct the mapping matrix G [Equation 13]
	G = interpG(N,xF,M,xC)
	GT = transpose(G)

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

	# Construct the coarse-grid iteration matrix [Equation 15]
	itMatCoarse = GT * itMatFull * G
	itMatCoarse = lu(itMatCoarse)

	# Pre-allocate full Nx1 solutions
	c = zeros(N,K+1)
	c[:,1] = par.c0(xF)

	# Pre-allocate coarse Mx1 solutions
	C = zeros(M,K+1)
	C[:,1] = par.c0(xC)

	# Pre-allocate RHS vectors for the full (stepTmpbFull) NxN problem, and the coarse (stepTmpbCoarse) MxM problem
	stepTmpbFull = similar(xF)
	stepTmpbCoarse = similar(xC)

	# Time stepping
	for k in 1:K

		# Construct the full Nx1 RHS vector
		stepTmpbFull .= c[:,k]
		stepTmpbFull .+= bval[:,k+1]

		# Construct the coarse Mx1 RHS vector
		stepTmpbCoarse .= GT * stepTmpbFull

		# Approximate solution on the coarse grid [Equation 15]
		ldiv!(itMatCoarse,stepTmpbCoarse)
		C[:,k+1] .= stepTmpbCoarse

		# Reconstruct solution on the fine grid [Equation 14]
		c[:,k+1] .= G * C[:,k+1]

	end

	return c, C

end
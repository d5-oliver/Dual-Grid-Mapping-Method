using LinearAlgebra

include("fvmStd.jl")

function FVM(par)

	# Transport coefficients and problem parameters
	x0 = par.x0
	xL = par.xL
	N = par.N
	xF = par.xF

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

	# Construct the iteration matrix A [Equation 7] and vector b [Equation 8]
	A, b = fvmStd(xF, R, D, v, mu, Gamma, omega, alpha0, beta0, sigma0, alphaL, betaL, sigmaL)
	itMat = I - tau * A
	itMat = lu(itMat)

	# Pre-allocate the b vectors for use in the time-stepping loop
	bval = zeros(N,K+1)
	bval[1,:] .= tau .* (b[1] .* g0 .+ Gamma(x0) / R(x0))
	bval[N,:] .= tau .* (b[N] .* gL .+ Gamma(xL) / R(xL))
	@views @inbounds for k in 1:K
		bval[2:N-1,k] .= tau .* b[2:N-1]
	end

	# Pre-allocate solutions
	c = zeros(N, K+1)
	c[:,1] .= par.c0(xF)
	
	# Pre-allocate the RHS vector
	stepTmpb = similar(xF)

	# Time stepping
	@views @inbounds for k in 1:K

		# Construct the RHS vector
		stepTmpb .= c[:,k]
		stepTmpb .+= bval[:,k+1]

		# Approximate solution on the fine grid [Equation 9]
		ldiv!(itMat,stepTmpb)
		c[:,k+1] .= stepTmpb

	end

	return c

end
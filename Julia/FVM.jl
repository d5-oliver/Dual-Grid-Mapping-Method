using LinearAlgebra

include("fvmStd.jl")

function FVM(par)

	x0::Int64 = par.x0
	xL::Int64 = par.xL
	N::Int64 = par.N
	xF::StepRangeLen{Float64} = par.xF

	tau::Float64 = par.tau
	t::StepRangeLen{Float64} = par.t
	K::Int64 = par.K

	R::Any = par.R
	D::Any = par.D
	v::Any = par.v
	mu::Any = par.mu
	Gamma::Any = par.Gamma

	omega::Int64 = par.omega

	alpha0::Float64 = par.alpha0
	beta0::Float64 = par.beta0
	sigma0::Float64 = par.sigma0
	g0::Any = par.g0(t)

	alphaL::Float64 = par.alphaL
	betaL::Float64 = par.betaL
	sigmaL::Float64 = par.sigmaL
	gL::Any = par.gL(t)

	A::Array{Float64}, b::Array{Float64} = fvmStd(xF,R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL)
	itMat::Array{Float64} = I - tau*A

	b0::Float64 = b[1]
	bL::Float64 = b[N]

	b[2:N-1] = tau*b[2:N-1]

	val0::Float64 = Gamma(x0)/R(x0)
	valL::Float64 = Gamma(xL)/R(xL)

	c::Array{Float64} = zeros(N,K+1)
	c[:,1] = par.c0(xF)

	for k in 1:K

		b[1] = tau * (b0*g0[k+1] + val0)
		b[N] = tau * (bL*gL[k+1] + valL)

		c[:,k+1] = itMat \ (c[:,k] + b)

	end

	return c

end
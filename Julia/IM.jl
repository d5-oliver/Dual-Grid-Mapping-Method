using LinearAlgebra

include("fvmStd.jl")

function interpG(N,xF,M,xC)

	G::Array{Float64} = zeros(N,M)
    P::Int64 = (N-1)/(M-1) + 1
    idx::StepRangeLen{Int64} = 0:0

    for m in 1:M-1

        idx = (m-1)*(P-1)+1:m*(P-1)+1

        G[idx,m] = (xC[m+1] .- xF[idx]) ./ (xC[m+1] .- xC[m])
        G[idx,m+1] = (xF[idx] .- xC[m]) ./ (xC[m+1] .- xC[m])

    end

	return G

end

function IM(par)

	x0::Float64 = par.x0
	xL::Float64 = par.xL

	N::Int64 = par.N
	M::Int64 = par.M

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

	alpha0::Float64 = par.alpha0
	beta0::Float64 = par.beta0
	sigma0::Float64 = par.sigma0
	g0::Any = par.g0(t)

	alphaL::Float64 = par.alphaL
	betaL::Float64 = par.betaL
	sigmaL::Float64 = par.sigmaL
	gL::Any = par.gL(t)

    G::Array{Float64} = interpG(N,xF,M,xC)
    GT::Array{Float64} = transpose(G)

	A::Array{Float64}, b::Array{Float64} = fvmStd(xF,R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL)
	itMatFull::Array{Float64} = I - tau*A

	b0::Float64 = b[1]
	bL::Float64 = b[N]

	b[2:N-1] = tau*b[2:N-1]

    val0::Float64 = Gamma(x0)/R(x0)
	valL::Float64 = Gamma(xL)/R(xL)

	itMatCoarse::Array{Float64} = GT * itMatFull * G

    c::Array{Float64} = zeros(N,K+1)
	c[:,1] = par.c0(xF)

	C::Array{Float64} = zeros(M,K+1)
	C[:,1] = par.c0(xC)

    for k in 1:K

        b[1] = tau*(b0*g0[k+1] + val0)
        b[N] = tau*(bL*gL[k+1] + valL)

        C[:,k+1] = itMatCoarse \ (GT * (c[:,k] + b))

        c[:,k+1] = G * C[:,k+1]

    end

	return c, C

end
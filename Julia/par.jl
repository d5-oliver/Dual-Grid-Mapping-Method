struct par

	M::Int64
	N::Int64

	x0::Float64
	xL::Float64

	xC::StepRangeLen{Float64}
	xF::StepRangeLen{Float64}

	T::Float64
	K::Int64
	tau::Float64
	t::StepRangeLen{Float64}

	R::Any
	D::Any
	v::Any
	mu::Any
	Gamma::Any

	omega::Int8

	sigma::Float64

	alpha0::Float64
	beta0::Float64
	sigma0::Float64
	g0::Any

	alphaL::Float64
	betaL::Float64
	sigmaL::Float64
	gL::Any

	c0::Any

	tSteps::Vector{Int64}

end

# Basic Heaviside step function.
function heaviside(x)

	return .!(x .< 0)

end

# Basic piecewise-constant function.
function layerEval(x::Any, layers::Array{Float64}, parcol::Int64)

	local y::Union{Int64,Float64,Array{Float64}}
	
	N::Int64 = length(x)

	if N == 1

		y = 0

	else

		y = zeros(N)

	end

	for ii in 1:size(layers,1)-1

		y = y .+ layers[ii,parcol].*heaviside(x .- layers[ii,1]) .- layers[ii,parcol].*heaviside(x .- layers[ii+1,1]);

	end

	return y

end

# Set test problem values.
# Cases 1 and 2 were defined from test problems used in:
# Carr, E.J. 2020. New semi-analytical solutions for advection-dispersion 
# equations in multilayer porous media. Transport in Porous Media, 
# 135(1):pp. 39–58. https://dx.doi.org/10.1007/s11242-020-01468-z.
function constructPar(x0::Float64, xL::Float64, T::Float64, N::Int64=1001, M::Int64=11, K::Int64=3600, omega::Int64=1, sigma::Float64=1e+10, testCase::Int64=1)

	local layers::Array{Float64}
	
	local xF::StepRangeLen{Float64}
	local xC::StepRangeLen{Float64}

	local tau::Float64

	local R::Any
	local D::Any
	local v::Any
	local mu::Any
	local Gamma::Any

	local alpha0::Float64
	local beta0::Float64
	local sigma0::Float64
	local g0::Any

	local alphaL::Float64
	local betaL::Float64
	local sigmaL::Float64
	local gL::Any

	local c0::Any

	local tSteps::Array{Int64}

	if (testCase == 1)

		layers = [
				-1.0 1.7 4.0 2.8 0.0 0.0
				10.0 7.0 4.0 9.0 0.0 0.0
				12.0 1.7 4.0 2.8 0.0 0.0
				20.0 7.0 4.0 9.0 0.0 0.0
				22.0 1.7 4.0 2.8 0.0 0.0
				31.0 1.7 4.0 2.8 0.0 0.0
			]

		xF = x0:xL/(N-1):xL
		xC = x0:xL/(M-1):xL

		tau = T/K
		t = 0:tau:T

		R = x -> layerEval(x,layers,2)
		D = x -> layerEval(x,layers,4)
		v = x -> layerEval(x,layers,3)
		mu = x -> layerEval(x,layers,5)
		Gamma = x -> layerEval(x,layers,6)

		alpha0 = v(x0) .- 1
		beta0 = D(x0)
		sigma0 = 1
		g0 = t -> v(x0)*ones(size(t))

		alphaL = v(xL)
		betaL = D(xL)
		sigmaL = v(xL)
		gL = t -> zeros(size(t))

		c0 = x -> zeros(size(x))

		tSteps = round.(K*[0.1,1,2,6,10,14,18]/T)

	elseif (testCase == 2)

		layers = [
			-1.0 1.7 4.0 2.8 0.0 0.0
			10.0 7.0 4.0 9.0 0.0 0.0
			12.0 1.7 4.0 2.8 0.0 0.0
			20.0 7.0 4.0 9.0 0.0 0.0
			22.0 1.7 4.0 2.8 0.0 0.0
			31.0 1.7 4.0 2.8 0.0 0.0
		]

		xF = x0:xL/(N-1):xL
		xC = x0:xL/(M-1):xL

		tau =  T/K
		t =  0:tau:T

		R = x -> layerEval(x,layers,2)
		D = x -> layerEval(x,layers,4)
		v = x -> layerEval(x,layers,3)
		mu = x -> layerEval(x,layers,5)
		Gamma = x -> layerEval(x,layers,6)

		alpha0 = v(x0) .- 1
		beta0 = D(x0)
		sigma0 = 1
		g0 = t -> v(x0)*heaviside(3 .- t)

		alphaL = v(xL)
		betaL = D(xL)
		sigmaL = v(xL)
		gL = t -> zeros(size(t))

		c0 = x -> zeros(size(x))

		tSteps = round.(K*[0.1,1,2,6,10,14,18]/T)

	elseif (testCase == 3)

		layers = [
			-1.0 1.7 4.0 2.8 0.0 0.0
			 0.5 7.0 4.0 9.0 0.0 0.0
			 1.0 1.7 4.0 2.8 0.0 0.0
			 1.5 7.0 4.0 9.0 0.0 0.0
			 2.0 1.7 4.0 2.8 0.0 0.0
			 2.5 7.0 4.0 9.0 0.0 0.0
			 3.0 1.7 4.0 2.8 0.0 0.0
			 3.5 7.0 4.0 9.0 0.0 0.0
			 4.0 1.7 4.0 2.8 0.0 0.0
			 4.5 7.0 4.0 9.0 0.0 0.0
			 5.0 1.7 4.0 2.8 0.0 0.0
			 5.5 7.0 4.0 9.0 0.0 0.0
			 6.0 1.7 4.0 2.8 0.0 0.0
			 6.5 7.0 4.0 9.0 0.0 0.0
			 7.0 1.7 4.0 2.8 0.0 0.0
			 7.5 7.0 4.0 9.0 0.0 0.0
			 8.0 1.7 4.0 2.8 0.0 0.0
			 8.5 7.0 4.0 9.0 0.0 0.0
			 9.0 1.7 4.0 2.8 0.0 0.0
			 9.5 7.0 4.0 9.0 0.0 0.0
			10.0 1.7 4.0 2.8 0.0 0.0
			10.5 7.0 4.0 9.0 0.0 0.0
			11.0 1.7 4.0 2.8 0.0 0.0
			11.5 7.0 4.0 9.0 0.0 0.0
			12.0 1.7 4.0 2.8 0.0 0.0
			12.5 7.0 4.0 9.0 0.0 0.0
			13.0 1.7 4.0 2.8 0.0 0.0
			13.5 7.0 4.0 9.0 0.0 0.0
			14.0 1.7 4.0 2.8 0.0 0.0
			14.5 7.0 4.0 9.0 0.0 0.0
			15.0 1.7 4.0 2.8 0.0 0.0
			15.5 7.0 4.0 9.0 0.0 0.0
			16.0 1.7 4.0 2.8 0.0 0.0
			16.5 7.0 4.0 9.0 0.0 0.0
			17.0 1.7 4.0 2.8 0.0 0.0
			17.5 7.0 4.0 9.0 0.0 0.0
			18.0 1.7 4.0 2.8 0.0 0.0
			18.5 7.0 4.0 9.0 0.0 0.0
			19.0 1.7 4.0 2.8 0.0 0.0
			19.5 7.0 4.0 9.0 0.0 0.0
			20.0 1.7 4.0 2.8 0.0 0.0
			20.5 7.0 4.0 9.0 0.0 0.0
			21.0 1.7 4.0 2.8 0.0 0.0
			21.5 7.0 4.0 9.0 0.0 0.0
			22.0 1.7 4.0 2.8 0.0 0.0
			22.5 7.0 4.0 9.0 0.0 0.0
			23.0 1.7 4.0 2.8 0.0 0.0
			23.5 7.0 4.0 9.0 0.0 0.0
			24.0 1.7 4.0 2.8 0.0 0.0
			24.5 7.0 4.0 9.0 0.0 0.0
			25.0 1.7 4.0 2.8 0.0 0.0
			25.5 7.0 4.0 9.0 0.0 0.0
			26.0 1.7 4.0 2.8 0.0 0.0
			26.5 7.0 4.0 9.0 0.0 0.0
			27.0 1.7 4.0 2.8 0.0 0.0
			27.5 7.0 4.0 9.0 0.0 0.0
			28.0 1.7 4.0 2.8 0.0 0.0
			28.5 7.0 4.0 9.0 0.0 0.0
			29.0 1.7 4.0 2.8 0.0 0.0
			29.5 7.0 4.0 9.0 0.0 0.0
			31.0 1.7 4.0 2.8 0.0 0.0
		]

		xF = x0:xL/(N-1):xL
		xC = x0:xL/(M-1):xL

		tau =  T/K
		t =  0:tau:T

		R = x -> layerEval(x,layers,2)
		D = x -> layerEval(x,layers,4)
		v = x -> layerEval(x,layers,3)
		mu = x -> layerEval(x,layers,5)
		Gamma = x -> layerEval(x,layers,6)

		alpha0 = v(x0) .- 1
		beta0 = D(x0)
		sigma0 = 1
		g0 = t -> v(x0)*ones(size(t))

		alphaL = v(xL)
		betaL = D(xL)
		sigmaL = v(xL)
		gL = t -> zeros(size(t))

		c0 = x -> zeros(size(x))

		tSteps = round.(K*[0.1,1,2,6,10,14,18]/T)

	elseif (testCase == 4)

		layers = [
			-1.0 1.7 4.0 2.8 0.0 0.0
			 0.5 7.0 4.0 9.0 0.0 0.0
			 1.0 1.7 4.0 2.8 0.0 0.0
			 1.5 7.0 4.0 9.0 0.0 0.0
			 2.0 1.7 4.0 2.8 0.0 0.0
			 2.5 7.0 4.0 9.0 0.0 0.0
			 3.0 1.7 4.0 2.8 0.0 0.0
			 3.5 7.0 4.0 9.0 0.0 0.0
			 4.0 1.7 4.0 2.8 0.0 0.0
			 4.5 7.0 4.0 9.0 0.0 0.0
			 5.0 1.7 4.0 2.8 0.0 0.0
			 5.5 7.0 4.0 9.0 0.0 0.0
			 6.0 1.7 4.0 2.8 0.0 0.0
			 6.5 7.0 4.0 9.0 0.0 0.0
			 7.0 1.7 4.0 2.8 0.0 0.0
			 7.5 7.0 4.0 9.0 0.0 0.0
			 8.0 1.7 4.0 2.8 0.0 0.0
			 8.5 7.0 4.0 9.0 0.0 0.0
			 9.0 1.7 4.0 2.8 0.0 0.0
			 9.5 7.0 4.0 9.0 0.0 0.0
			10.0 1.7 4.0 2.8 0.0 0.0
			10.5 7.0 4.0 9.0 0.0 0.0
			11.0 1.7 4.0 2.8 0.0 0.0
			11.5 7.0 4.0 9.0 0.0 0.0
			12.0 1.7 4.0 2.8 0.0 0.0
			12.5 7.0 4.0 9.0 0.0 0.0
			13.0 1.7 4.0 2.8 0.0 0.0
			13.5 7.0 4.0 9.0 0.0 0.0
			14.0 1.7 4.0 2.8 0.0 0.0
			14.5 7.0 4.0 9.0 0.0 0.0
			15.0 1.7 4.0 2.8 0.0 0.0
			15.5 7.0 4.0 9.0 0.0 0.0
			16.0 1.7 4.0 2.8 0.0 0.0
			16.5 7.0 4.0 9.0 0.0 0.0
			17.0 1.7 4.0 2.8 0.0 0.0
			17.5 7.0 4.0 9.0 0.0 0.0
			18.0 1.7 4.0 2.8 0.0 0.0
			18.5 7.0 4.0 9.0 0.0 0.0
			19.0 1.7 4.0 2.8 0.0 0.0
			19.5 7.0 4.0 9.0 0.0 0.0
			20.0 1.7 4.0 2.8 0.0 0.0
			20.5 7.0 4.0 9.0 0.0 0.0
			21.0 1.7 4.0 2.8 0.0 0.0
			21.5 7.0 4.0 9.0 0.0 0.0
			22.0 1.7 4.0 2.8 0.0 0.0
			22.5 7.0 4.0 9.0 0.0 0.0
			23.0 1.7 4.0 2.8 0.0 0.0
			23.5 7.0 4.0 9.0 0.0 0.0
			24.0 1.7 4.0 2.8 0.0 0.0
			24.5 7.0 4.0 9.0 0.0 0.0
			25.0 1.7 4.0 2.8 0.0 0.0
			25.5 7.0 4.0 9.0 0.0 0.0
			26.0 1.7 4.0 2.8 0.0 0.0
			26.5 7.0 4.0 9.0 0.0 0.0
			27.0 1.7 4.0 2.8 0.0 0.0
			27.5 7.0 4.0 9.0 0.0 0.0
			28.0 1.7 4.0 2.8 0.0 0.0
			28.5 7.0 4.0 9.0 0.0 0.0
			29.0 1.7 4.0 2.8 0.0 0.0
			29.5 7.0 4.0 9.0 0.0 0.0
			31.0 1.7 4.0 2.8 0.0 0.0
		]

		xF = x0:xL/(N-1):xL
		xC = x0:xL/(M-1):xL

		tau =  T/K
		t =  0:tau:T

		R = x -> layerEval(x,layers,2)
		D = x -> layerEval(x,layers,4)
		v = x -> layerEval(x,layers,3)
		mu = x -> layerEval(x,layers,5)
		Gamma = x -> layerEval(x,layers,6)

		alpha0 = v(x0) .- 1
		beta0 = D(x0)
		sigma0 = 1
		g0 = t -> v(x0)*heaviside(3 .- t)

		alphaL = v(xL)
		betaL = D(xL)
		sigmaL = v(xL)
		gL = t -> zeros(size(t))

		c0 = x -> zeros(size(x))

		tSteps = round.(K*[0.1,1,2,6,10,14,18]/T)

	elseif (testCase == 5)

		xF = x0:xL/(N-1):xL
		xC = x0:xL/(M-1):xL

		tau =  T/K
		t =  0:tau:T

		R = x -> 0.75 .+ 0.25*sin.(5*x)
		D = x -> 14.85 .+ 12.55*cos.(5*x)
		v = x -> 10 .+ zeros(size(x))
		mu = x -> 0.05 .+ 0.02*sin.(5*x)
		Gamma = x -> 0.025 .+ 0.01*sin.(5*x)

		alpha0 = v(x0)
		beta0 = D(x0)
		sigma0 = 1e+10
		g0 = t -> 1 .+ zeros(size(t))

		alphaL = v(xL)
		betaL = D(xL)
		sigmaL = v(xL)
		gL = t -> zeros(size(t))

		c0 = x -> zeros(size(x))

		tSteps = round.(K*([0.1,1,2,6,10,14,18]/10)/T)

	elseif (testCase == 6)

		xF = x0:xL/(N-1):xL
		xC = x0:xL/(M-1):xL

		tau =  T/K
		t =  0:tau:T

		local aheavi::Any = (k,x) -> 1 ./ (1 .+ exp.((-2)*k*x))
		R = x -> 0.75 .+ (1 .+ 0.75*sin.(4*x)) .* aheavi(2,x .- 10)
		D = x -> 9 .+ 0.3 * x .* cos.(3*x)
		v = x -> 14 .+ zeros(size(x))
		mu = x -> (0.25 .+ 0.15*sin.(2.5*x)) .* (1 .- aheavi(2,x .- 10))
		Gamma = x -> 0.01*(1 .+ sin.(x)) .* aheavi(2,x .- 10)

		alpha0 = v(x0) .- 1
		beta0 = D(x0)
		sigma0 = 1
		g0 = t -> v(x0)*heaviside(3 .- t)

		alphaL = v(xL)
		betaL = D(xL)
		sigmaL = v(xL)
		gL = t -> zeros(size(t))

		c0 = x -> zeros(size(x))

		tSteps = round.(K*[0.1,0.4,0.9,2,4,5,6]/T)

	end

	return par(M,N,x0,xL,xC,xF,T,K,tau,t,R,D,v,mu,Gamma,omega,sigma,alpha0,beta0,sigma0,g0,alphaL,betaL,sigmaL,gL,c0,tSteps)

end
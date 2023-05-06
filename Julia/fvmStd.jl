function fvmStd(x,R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL)

	N = length(x)

	dx = diff(x)
	Dx = similar(x)
	cvEW = similar(x,N+1)

	Di = similar(cvEW)
	mui = mu(x)
	Ri = R(x)
	Gammai = Gamma(x)

	Dx[1] = dx[1]/2
	for ii in 2:N-1
		Dx[ii] = (dx[ii-1] + dx[ii])/2
	end
	Dx[N] = dx[N-1]/2

	rho = @. 1 / (Ri * Dx)

	cvEW[1] = x[1]
	Di[1] = D(cvEW[1])
	for ii in 2:N
		cvEW[ii] = (x[ii] + x[ii-1])/2
		Di[ii] = D(cvEW[ii])/dx[ii-1]
	end
	cvEW[N+1] = x[N]
	Di[N+1] = D(cvEW[N+1])

	vi = v(cvEW)

	A = zeros(N,N)

	b = similar(x)

	# Boundary node x=0
	A[1,1] = -rho[1] * (Di[2] + vi[2] * (1 - 0.5*omega) + Di[1] * (sigma0 + alpha0)/beta0 - vi[1]) - mui[1]/Ri[1]
	A[1,2] = rho[1] * (Di[2] - vi[2] * 0.5*omega)

	b[1] = rho[1] * Di[1] * sigma0 / beta0

	# Interior nodes
	for ii in 2:N-1

		A[ii,ii-1] = rho[ii] * (Di[ii] + vi[ii] * (1 - 0.5*omega))
		A[ii,ii] = -rho[ii] * (Di[ii+1] + vi[ii+1] * (1 - 0.5*omega) + Di[ii] - vi[ii] * 0.5*omega) - mui[ii]/Ri[ii]
		A[ii,ii+1] = rho[ii] * (Di[ii+1] - vi[ii+1] * 0.5*omega)

		b[ii] = Gammai[ii] / Ri[ii]

	end

	# Boundary node x=L
	A[N,N-1] = rho[N] * (Di[N] + vi[N] * (1 - 0.5*omega))
	A[N,N] = -rho[N] * (Di[N+1] * (sigmaL - alphaL)/betaL + vi[N+1] + Di[N] - vi[N] * 0.5*omega) - mui[N]/Ri[N]

	b[N] = rho[N] * Di[N+1] * sigmaL / betaL

	return A, b

end
function constructArrays(R,D,v,mu,Gamma,omega,sigma,xF,xC)

	N = length(xF)
	M = length(xC)
	P = div(N-1,M-1) + 1

	vi = v(xF)
	mui = mu(xF)
	Ri = R(xF)
	Gammai = Gamma(xF)

	A = zeros(N+M-2,N+M-2)
	b = zeros(P,M-1)
	uv = zeros(P,2,M-1)

	# Concatenate local iteration matrices in a single block diagonal matrix
	@views @inbounds for m in 1:M-1

		idxGrid = (m-1)*(P-1)+1:m*(P-1)+1
		idxStore = (m-1)*P+1:m*P

		fvmSig!(A[idxStore,idxStore],b[:,m],uv[:,:,m],P,xF[idxGrid],Ri[idxGrid],D,vi[idxGrid],mui[idxGrid],Gammai[idxGrid],omega,sigma)

	end

	return A, b, uv

end

function fvmSig!(Am,bm,uvm,P,xF,Ri,D,vi,mui,Gammai,omega,sigma)

	dx = diff(xF)

	Dx = similar(xF)
	Dx[1] = dx[1]/2
	for ii in 2:P-1
		Dx[ii] = (dx[ii-1] + dx[ii])/2
	end
	Dx[P] = dx[P-1]/2

	cvEW = similar(xF,P-1)
	Di = similar(xF,P-1)
	for ii in 1:P-1
		cvEW[ii] = (xF[ii+1] + xF[ii])/2
		Di[ii] = D(cvEW[ii])/dx[ii]
	end

	rho = 1 ./ (Ri .* Dx)

	bm .= Gammai ./ Ri

	uvm[1,1] = rho[1] * sigma
	uvm[P,2] = rho[P] * sigma

	# Boundary node x=Xm
	Am[1,1] = -rho[1] * (Di[1] + vi[1] * (1 - omega/2) + sigma) - mui[1]/Ri[1]
	Am[1,2] = rho[1] * (Di[1] - vi[1] * omega/2)

	# Interior nodes
	@inbounds for ii in 2:P-1

		Am[ii,ii-1] = rho[ii] * (Di[ii-1] + vi[ii-1] * (1 - omega/2))
		Am[ii,ii] = -rho[ii] * (Di[ii] + vi[ii] * (1 - omega/2) + Di[ii-1] - vi[ii-1] * omega/2) - mui[ii]/Ri[ii]
		Am[ii,ii+1] = rho[ii] * (Di[ii] - vi[ii] * (1 - omega/2))

	end

	# Boundary node x=Xm+1
	Am[P,P-1] = rho[P] * (Di[P-1] + vi[P-1] * (1 - omega/2))
	Am[P,P] = -rho[P] * (Di[P-1] - vi[P-1] * omega/2 + sigma) - mui[P]/Ri[P]

end

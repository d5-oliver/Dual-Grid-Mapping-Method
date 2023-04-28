function constructArrays(R,D,v,mu,Gamma,omega,sigma,xF,xC)

	N::Int64 = length(xF)
	M::Int64 = length(xC)
	P::Int64 = (N-1)/(M-1) + 1

	vi::Array{Float64} = v(xF)
	mui::Array{Float64} = mu(xF)
	Ri::Array{Float64} = R(xF)
	Gammai::Array{Float64} = Gamma(xF)

	A::Array{Float64} = zeros(N+M-2,N+M-2)
	b::Array{Float64} = zeros(P,M-1)
	d::Array{Float64} = zeros(P+1,M-1)

	idxGrid::StepRangeLen{Int64} = 1:P
	idxStore::StepRangeLen{Int64} = 1:P

	# Concatenate local iteration matrices in a single block diagonal matrix
	for m in 1:M-1

		idxGrid = (m-1)*(P-1)+1:m*(P-1)+1;
		idxStore = (m-1)*P+1:m*P;

		A[idxStore,idxStore], b[:,m], d[:,m] = fvmSig(P,xF[idxGrid],Ri[idxGrid],D,vi[idxGrid],mui[idxGrid],Gammai[idxGrid],omega,sigma)

	end

	return A, b, d

end

function fvmSig(P,xF,Ri,D,vi,mui,Gammai,omega,sigma)

	dx::Array{Float64} = diff(xF)

	Dx::Array{Float64} = zeros(P)
	Dx[1] = dx[1]/2
	for ii in 2:P-1
		Dx[ii] = (dx[ii-1] + dx[ii])/2
	end
	Dx[P] = dx[P-1]/2

	cvEW::Array{Float64} = zeros(P-1)
	Di::Array{Float64} = zeros(P-1)
	for ii in 1:P-1
		cvEW[ii] = (xF[ii+1] + xF[ii])/2
		Di[ii] = D(cvEW[ii])/dx[ii]
	end

	rho::Array{Float64} = 1 ./ (Ri .* Dx)

	bm::Array{Float64} = Gammai ./ Ri
	dm::Array{Float64} = zeros(P+1)
	dm[1] = rho[1] * sigma
	dm[P+1] = rho[P] * sigma
	Am::Array{Float64} = zeros(P,P)

	# Boundary node x=Xm
	Am[1,1] = -rho[1] * (Di[1] + vi[1] * (1 - omega/2) + sigma) - mui[1]/Ri[1]
	Am[1,2] = rho[1] * (Di[1] - vi[1] * omega/2)

	# Interior nodes
	for ii in 2:P-1

		Am[ii,ii-1] = rho[ii] * (Di[ii-1] + vi[ii-1] * (1 - omega/2))
		Am[ii,ii] = -rho[ii] * (Di[ii] + vi[ii] * (1 - omega/2) + Di[ii-1] - vi[ii-1] * omega/2) - mui[ii]/Ri[ii]
		Am[ii,ii+1] = rho[ii] * (Di[ii] - vi[ii] * (1 - omega/2))

	end

	# Boundary node x=Xm+1
	Am[P,P-1] = rho[P] * (Di[P-1] + vi[P-1] * (1 - omega/2))
	Am[P,P] = -rho[P] * (Di[P-1] - vi[P-1] * omega/2 + sigma) - mui[P]/Ri[P]

	return Am, bm, dm

end


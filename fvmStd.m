function [A,b] = fvmStd(R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL,x)
%FVMSTD constructs arrays A and b using the finite volume method.
%
% FVMSTD(R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL,x)
% constructs arrays for a standard vertex-centered finite volume method.
%
% Transport coefficients R, D, v, mu, and Gamma are function handles.
%
% Boundary coefficients alpha0, beta0, and sigma0 correspond to the
% boundary condition enforced on the left side of the domain.
%
% Boundary coefficients alphaL, betaL, and sigmaL correspond to the
% boundary condition enforced on the right side of the domain.
%
% Approximation parameter omega implements either averaging or upwinding in
% the finite volume scheme.

N = length(x);

dx = diff(x);
Dx = 0.5*[dx(1);dx(1:N-2) + dx(2:N-1);dx(N-1)];
cvEW = 0.5*(x(2:N) + x(1:N-1));

Di = D(cvEW) ./ dx;
vi = v(cvEW);
Ri = R(x);
mui = mu(x);
Gammai = Gamma(x);

rho = 1 ./ (Ri .* Dx);

A = zeros(N,N);


% Left boundary node
A(1,1) = -rho(1) * ( Di(1) + vi(1) * (1 - 0.5*omega) + D(x(1)) * (sigma0 + alpha0)/beta0 - v(x(1)) ) - mui(1)/Ri(1);
A(1,2) = rho(1) * ( Di(1) - vi(1) * 0.5*omega );

b0 = rho(1) * D(x(1)) * sigma0 / beta0;


% Interior nodes
A(2:N+1:N*(N-2)-1) = rho(2:N-1) .* (Di(1:N-2) + vi(1:N-2) * (1 - 0.5*omega));
A(N+2:N+1:N*(N-1)-1) = -rho(2:N-1) .* ( Di(2:N-1) + vi(2:N-1) * (1 - 0.5*omega) + Di(1:N-2) - vi(1:N-2) * 0.5*omega ) - mui(2:N-1)./Ri(2:N-1);
A(2*N+2:N+1:N*N-1) = rho(2:N-1) .* (Di(2:N-1) - vi(2:N-1) * 0.5*omega);


% Right boundary node
A(N,N-1) = rho(N) * (Di(N-1) + vi(N-1) * (1 - 0.5*omega));
A(N,N) = -rho(N) * ( D(x(N))*(sigmaL - alphaL)/betaL + v(x(N)) + Di(N-1) - vi(N-1) * 0.5*omega ) - mui(N)/Ri(N);

bL = rho(N) * D(x(N)) * sigmaL / betaL;


% Define b vector
b = [b0;Gammai(2:N-1)./Ri(2:N-1);bL];

end

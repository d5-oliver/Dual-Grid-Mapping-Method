function [A,b,d] = fvmSig(R,D,v,mu,Gamma,omega,sigma,x,X)
%FVMSIG constructs all local arrays Am and bm using the finite volume 
%method.
%
% FVMSIG(R,D,v,mu,Gamma,omega,sigma,x,X) constructs arrays for a
% vertex-centered finite volume method.
%
% Transport coefficients R, D, v, mu, and Gamma are function handles.
%
% Boundaries utilise approximate Dirichlet conditions where the constant 
% sigma is very large.
%
% Approximation parameter omega implements either averaging or upwinding in
% the finite volume scheme.

N = length(x);
M = length(X);
P = (N-1) / (M-1) + 1;

A = zeros(N,N);
b = zeros(P,M-1);
d = zeros(P+1,M-1);

% Utilise a vertex-centered finite volume method to construct local arrays
% for the problem defined in [Equations 16-19].
for jj = 1:M-1

    idxGrid = (jj-1)*(P-1)+1:jj*(P-1)+1;
    idxStore = (jj-1)*P+1:jj*P;

    xLocal = x(idxGrid);

    dxLocal = diff(xLocal);
    DxLocal = 0.5*[dxLocal(1);dxLocal(1:P-2) + dxLocal(2:P-1);dxLocal(P-1)];
    cvEW = 0.5*(xLocal(2:P) + xLocal(1:P-1));

    Di = D(cvEW) ./ dxLocal;
    vi = v(cvEW);
    Ri = R(xLocal);
    mui = mu(xLocal);
    Gammai = Gamma(xLocal);

    rho = 1 ./ (Ri .* DxLocal);

    % Left boundary node
    A(idxStore(1),idxStore(1)) = -rho(1) * ( Di(1) + vi(1) * (1 - 0.5*omega) + sigma ) - mui(1)/Ri(1);
    A(idxStore(1),idxStore(2)) = rho(1) * ( Di(1) - vi(1) * 0.5*omega );

    d(1,jj) = rho(1) * sigma;

    % Interior nodes
    for ii = 2:P-1

        A(idxStore(ii),idxStore(ii-1)) = rho(ii) .* (Di(ii-1) + vi(ii-1) * (1 - 0.5*omega));
        A(idxStore(ii),idxStore(ii)) = -rho(ii) .* ( Di(ii) + vi(ii) * (1 - 0.5*omega) + Di(ii-1) - vi(ii-1) * 0.5*omega ) - mui(ii)./Ri(ii);
        A(idxStore(ii),idxStore(ii+1)) = rho(ii) .* (Di(ii) - vi(ii) * 0.5*omega);

    end

    % Right boundary node
    A(idxStore(P),idxStore(P-1)) = rho(P) * ( Di(P-1) + vi(P-1) * (1 - 0.5*omega) );
    A(idxStore(P),idxStore(P)) = -rho(P) * ( Di(P-1) - vi(P-1) * 0.5*omega + sigma ) - mui(P)/Ri(P);

    d(end,jj) = rho(P) * sigma;

    % Local b vector
    b(:,jj) = Gammai ./ Ri;

end

end

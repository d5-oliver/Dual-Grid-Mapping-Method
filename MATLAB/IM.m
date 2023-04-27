function [C,c] = IM(par)
%IM utilises an interpolation-based mapping to solve the 1D ADR equation.
%
% IM(par) assumes par is a struct array containing information about the 
% transport problem. 

% Transport coefficients and problem parameters
N = par.N;

M = par.M;

xF = par.xF(:);
xC = par.xC(:);

t0 = par.t0;
tK = par.tK;
K = par.K;
tau = (tK - t0) / K;
t = t0:tau:tK;

R = par.R;
D = par.D;
v = par.v;
mu = par.mu;
Gamma = par.Gamma;

omega = par.omega;

alpha0 = par.alpha0;
beta0 = par.beta0;
sigma0 = par.sigma0;

alphaL = par.alphaL;
betaL = par.betaL;
sigmaL = par.sigmaL;

g0 = par.g0(t);
gL = par.gL(t);

C = zeros(M,K+1);
C(:,1) = par.c0(xC);

c = zeros(N,K+1);
c(:,1) = par.c0(xF);

I = eye(N,N);

% Construct the mapping matrix G [Equation 13]
G = interpG(xF,xC);

% Construct the iteration matrix A [Equation 7] and vector b [Equation 8]
[A,b] = fvmStd(R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL,xF);

b0 = b(1);
bL = b(N);

b(2:N-1) = tau*b(2:N-1);

val0 = Gamma(xF(1)) / R(xF(1));
valL = Gamma(xF(N)) / R(xF(N));

itMat = G' * (I - tau*A) * G;

% Time stepping
for kk = 1:K

    b(1) = tau * (b0*g0(kk+1) + val0);

    b(N) = tau * (bL*gL(kk+1) + valL);

    % Approximate solution on the coarse grid [Equation 15]
    C(:,kk+1) = itMat \ ( G' * ( c(:,kk) + b ) );

    % Reconstruct solution on the fine grid [Equation 14]
    c(:,kk+1) = G * C(:,kk+1);

end

end

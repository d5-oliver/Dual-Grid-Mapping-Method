function c = FVM(par)
%FVM utilises the finite volume method to solve the 1D ADR equation.
%
% FVM(par) assumes par is a struct array containing information about the
% transport problem. 

% Transport coefficients and problem parameters
N = par.N;

xF = par.xF(:);

tK = par.tK;
K = par.K;
tau = tK / K;
t = 0:tau:tK;

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

c = zeros(N,K+1);
c(:,1) = par.c0(xF);

I = eye(N,N);

% Construct the iteration matrix A [Equation 7] and vector b [Equation 8]
[A,b] = fvmStd(R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL,xF);

itMat = I - tau*A;

b0 = b(1);
bL = b(N);

b(2:N-1) = tau*b(2:N-1);

val0 = Gamma(xF(1)) / R(xF(1));
valL = Gamma(xF(N)) / R(xF(N));

% Time stepping
for k = 1:K

    b(1) = tau * (b0*g0(k+1) + val0);

    b(N) = tau * (bL*gL(k+1) + valL);
    
    % Approximate solution on the fine grid [Equation 9]
    c(:,k+1) = itMat \ ( c(:,k) + b );

end

end

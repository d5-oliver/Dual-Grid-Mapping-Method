function [C,c] = DM(par)
%DM utilises a discretisation-based mapping to solve the 1D ADR equation.
%
% DM(par) assumes par is a struct array containing information about the 
% transport problem. 

% Transport coefficients and problem parameters
N = par.N;

M = par.M;

xF = par.xF(:);
xC = par.xC(:);

tK = par.tK;
K = par.K;
tau = tK / K;
t = 0:tau:tK;

D = par.D;
v = par.v;
mu = par.mu;
Gamma = par.Gamma;

omega = par.omega;

sigma = par.sigma;

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

P = (N - 1) / (M - 1) + 1;

I = eye(N+M-2,N+M-2);

Gm = zeros(N+M-2,M);

phim = zeros(N+M-2,1);

idx = setdiff( 1:(N+M-2) , (1:M-2)*P+1 );

% Construct the iteration matrix A [Equation 7] and vector b [Equation 8]
[A,b] = fvmStd(par.R,D,v,mu,Gamma,omega,alpha0,beta0,sigma0,alphaL,betaL,sigmaL,xF);
itMatFull = I(1:N,1:N) - tau*A;

% Construct the iteration matrices A [Equation 21] and vectors 
% b, u, and v [Equation 23] corresponding to locally defined sub-problems 
% [Equations 16-19]
[Am,bm,uvm] = fvmSig(par.R,D,v,mu,Gamma,omega,sigma,xF,xC);
itMatm = I - tau*Am;

b0 = b(1);
bL = b(N);

b(2:N-1) = tau*b(2:N-1);

val0 = Gamma(xF(1)) / par.R(xF(1));
valL = Gamma(xF(N)) / par.R(xF(N));

% Construct the mapping matrix G [Equation 25]
for m = 1:M-1

    idxStore = (m-1)*P+1:m*P;

    Gm(idxStore,[m,m+1]) = itMatm(idxStore,idxStore) \ ( tau * [uvm(1:P,m),uvm(2:P+1,m)] );

end

itMatCoarse = Gm(idx,:)' * itMatFull * Gm(idx,:);

% Time stepping
for k = 1:K
    
    % Construct the phi vector [Equation 26]
    for m = 1:M-1

        idxGrid = (m-1)*(P - 1)+1:m*(P - 1)+1;

        idxStore = (m-1)*P+1:m*P;

        phim(idxStore) = itMatm(idxStore,idxStore) \ ( tau*bm(:,m) + c(idxGrid,k) );

    end

    b(1) = tau * (b0*g0(k+1) + val0);

    b(N) = tau * (bL*gL(k+1) + valL);

    % Approximate solution on the coarse grid [Equation 12]
    C(:,k+1) = itMatCoarse \ ( Gm(idx,:)' * ( c(:,k) - itMatFull * phim(idx) + b ) );

    % Reconstruct solution on the fine grid [Equation 11]
    c(:,k+1) = Gm(idx,:) * C(:,k+1) + phim(idx);

end

end

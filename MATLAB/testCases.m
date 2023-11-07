%% Define test case parameters
% Cases 1 and 2 were defined from test problems used in:
% Carr, E.J. 2020. New semi-analytical solutions for advection-dispersion
% equations in multilayer porous media. Transport in Porous Media,
% 135(1):pp. 39–58. https://dx.doi.org/10.1007/s11242-020-01468-z.

switch testCase

    % Case 9 from (Carr 2020)
    case 1

        % layers_testCase1and2.mat contains x-locations of interfaces
        % between material layers in the heterogeneous domain, and also
        % transport coefficient values for each corresponding layer.
        load('layers_testCase1and2.mat')

        % Transport coefficients
        par.R = @(x) layerEval(layers,2,x);
        par.v = @(x) layerEval(layers,3,x);
        par.D = @(x) layerEval(layers,4,x);
        par.mu = @(x) layerEval(layers,5,x);
        par.Gamma = @(x) layerEval(layers,6,x);

        % Boundary coefficients
        par.alpha0 = par.v(0) - 1;
        par.beta0 = par.D(par.xL);
        par.sigma0 = 1;

        par.alphaL = par.v(par.xL);
        par.betaL = par.D(par.xL);
        par.sigmaL = par.v(par.xL);

        par.g0 = @(t) par.v(0)*ones(size(t));
        par.gL = @(t) zeros(size(t));

        % Initial condition
        par.c0 = @(x) zeros(size(x));

        % Choose which values of 't' should be plotted
        tSteps = round(par.K*[0.1,1,2,6,10,14,18]/par.tK) + 1;

        % Case 10 from (Carr 2020)
    case 2

        load('layers_testCase1and2.mat')

        par.R = @(x) layerEval(layers,2,x);
        par.v = @(x) layerEval(layers,3,x);
        par.D = @(x) layerEval(layers,4,x);
        par.mu = @(x) layerEval(layers,5,x);
        par.Gamma = @(x) layerEval(layers,6,x);

        par.alpha0 = par.v(0) - 1;
        par.beta0 = par.D(0);
        par.sigma0 = 1;

        par.alphaL = par.v(par.xL);
        par.betaL = par.D(par.xL);
        par.sigmaL = par.v(par.xL);

        par.g0 = @(t) par.v(0)*heaviside(3-t);
        par.gL = @(t) zeros(size(t));

        par.c0 = @(x) zeros(size(x));

        tSteps = round(par.K*[0.1,1,2,6,10,14,18]/par.tK) + 1;

        % Test problem with many layers (similar to test problem 1)
    case 3

        load('layers_testCase3and4.mat')

        par.R = @(x) layerEval(layers,2,x);
        par.v = @(x) layerEval(layers,3,x);
        par.D = @(x) layerEval(layers,4,x);
        par.mu = @(x) layerEval(layers,5,x);
        par.Gamma = @(x) layerEval(layers,6,x);

        par.alpha0 = par.v(0) - 1;
        par.beta0 = par.D(0);
        par.sigma0 = 1;

        par.alphaL = par.v(par.xL);
        par.betaL = par.D(par.xL);
        par.sigmaL = par.v(par.xL);

        par.g0 = @(t) par.v(0)*ones(size(t));
        par.gL = @(t) zeros(size(t));

        par.c0 = @(x) zeros(size(x));

        tSteps = round(par.K*[0.1,1,2,6,10,14,18]/par.tK) + 1;

        % Test problem with many layers (similar to test problem 2)
    case 4

        load('layers_testCase3and4.mat')

        par.R = @(x) layerEval(layers,2,x);
        par.v = @(x) layerEval(layers,3,x);
        par.D = @(x) layerEval(layers,4,x);
        par.mu = @(x) layerEval(layers,5,x);
        par.Gamma = @(x) layerEval(layers,6,x);

        par.alpha0 = par.v(0) - 1;
        par.beta0 = par.D(0);
        par.sigma0 = 1;

        par.alphaL = par.v(par.xL);
        par.betaL = par.D(par.xL);
        par.sigmaL = par.v(par.xL);

        par.g0 = @(t) par.v(0)*heaviside(3-t);
        par.gL = @(t) zeros(size(t));

        par.c0 = @(x) zeros(size(x));

        tSteps = round(par.K*[0.1,1,2,6,10,14,18]/par.tK) + 1;

        % General new test problem - Similar to test problem 2
    case 5

        par.R = @(x) 0.75 + 0.25*sin(5*x);
        par.v = @(x) 10*ones(size(x));
        par.D = @(x) 14.85 + 12.55*cos(5*x);
        par.mu = @(x) 0.05 + 0.02*sin(5*x);
        par.Gamma = @(x) 0.025 + 0.01*sin(5*x);

        par.alpha0 = par.v(0);
        par.beta0 = par.D(0);
        par.sigma0 = 1e+10;

        par.alphaL = par.v(par.xL);
        par.betaL = par.D(par.xL);
        par.sigmaL = par.v(par.xL);

        par.g0 = @(t) ones(size(t));
        par.gL = @(t) zeros(size(t));

        par.c0 = @(x) zeros(size(x));

        tSteps = round(par.K*([0.1,1,2,6,10,14,18]/10)/par.tK) + 1;

        % Alternative general new test problem
    case 6

        aheavi = @(k,x) 1./(1 + exp(-2*k*x));
        par.R = @(x) 0.75 + (1 + 0.75*sin(4*x)) .* aheavi(2,x-10);
        par.v = @(x) 14*ones(size(x));
        par.D = @(x) 9 + 0.3 * x .* cos(3*x);
        par.mu = @(x) (0.25 + 0.15*sin(2.5*x)) .* (1 - aheavi(2,x-10));
        par.Gamma = @(x) 0.01*(1 + sin(x)) .* aheavi(2,x-10);

        par.alpha0 = par.v(0) - 1;
        par.beta0 = par.D(0);
        par.sigma0 = 1;

        par.alphaL = par.v(par.xL);
        par.betaL = par.D(par.xL);
        par.sigmaL = par.v(par.xL);

        par.g0 = @(t) par.v(0)*heaviside(3-t);
        par.gL = @(t) zeros(size(t));

        par.c0 = @(x) zeros(size(x));

        tSteps = round(par.K*[0.1,0.4,0.9,2,4,5,6]/par.tK) + 1;

        % Test problem with coefficients of varying magnitude by layer -
        % Similar to test problem 1
    case 7

        load('layers_testCase7.mat')

        % Transport coefficients
        par.R = @(x) layerEval(layers,2,x);
        par.v = @(x) layerEval(layers,3,x);
        par.D = @(x) layerEval(layers,4,x);
        par.mu = @(x) layerEval(layers,5,x);
        par.Gamma = @(x) layerEval(layers,6,x);

        % Boundary coefficients
        par.alpha0 = par.v(0) - 1;
        par.beta0 = par.D(0);
        par.sigma0 = 1;

        par.alphaL = par.v(par.xL);
        par.betaL = par.D(par.xL);
        par.sigmaL = par.v(par.xL);

        par.g0 = @(t) par.v(0)*ones(size(t));
        par.gL = @(t) zeros(size(t));

        % Initial condition
        par.c0 = @(x) zeros(size(x));

        % Choose which values of 't' should be plotted
        tSteps = round(par.K*[0.005,0.05,0.1,0.12,0.15,0.2,0.25]/par.tK) + 1;

end

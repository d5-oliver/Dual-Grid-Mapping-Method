%% Generate solutions for test cases from paper

%% Clear workspace, figures, and command window
close all
clearvars
clc

%% Choose solution and domain parameters
% Select test case (choices match figures, testCase = 1, 2, ..., 6)
testCase = 7;

% Choose number of spatial fine and coarse grid nodes.
par.N = 3601;
par.M = 11;

% Set domain, 0 <= x <= xL
par.xL = 30;

% Choose time step size. Size of a given time step is tK/K.
par.t0 = 0;
par.tK = 18;
par.K = 3600;

% Generate parameter set for chosen test case
run('testCases.m')

% If test case is a layered problem, ensure a node is placed on interfaces
if ismember(testCase,[1,2,3,4])
    while any(mod(layers(:,1),par.xL/(par.N - 1)))
        par.N = par.N + 1;
    end
end

% Define fine grid (xF)
par.xF = linspace(0,par.xL,par.N);

% Ensure coarse grid is a subset of fine grid nodes.
while mod( (par.N - 1)/(par.M - 1) , 1 ) > 0
    par.M = par.M + 1;
end

% Define coarse grid (xC)
par.xC = linspace(0,par.xL,par.M);

% Choose approximate-Dirichlet boundary condition scalar (see Equations 21
% and 22 in paper)
par.sigma = 1e+10;

% Select averaging or upstream weighting for approximate concentration at
% control volume faces.
% omega = 1 corresponds to averaging.
% omega = 0 corresponds to upwinding.
par.omega = 1;

% Choose to generate/save figures
plotFigs = false;
saveFigs = false;

%% Evaluation of different solution methods

tic
c_FVM = FVM(par); % Standard FVM solution
timerFVM = toc;

tic
[C_IM,c_IM] = IM(par); % Interpolation solution
timerIM = toc;

tic
[C_DM,c_DM] = DM(par); % Adaptive dual-grid solution
timerDM = toc;

disp(['Fine-grid solution (FVM) timer: ',num2str(timerFVM)])
disp(['Dual-grid mapping method (IM) timer: ',num2str(timerIM)])
disp(['Dual-grid mapping method (DM) timer: ',num2str(timerDM)])

%% Generate plots

if plotFigs

    cols = [
        0.6235    0.6549    0.6784
        0.1373    0.1529    0.1608
        0.0863    0.4431    0.8510
        0.7098    0.0353    0.0471
        ];

    figure;
    hold on
    box on

    % Fine grid solution plots, from both FVM and dual-grid (DM) methods
    for kk = tSteps
        plot(par.xF,c_FVM(:,kk),'Color',cols(1,:),'LineWidth',10)
        plot(par.xF,c_DM(:,kk),'--','Color',cols(2,:),'LineWidth',10)
    end

    % Coarse-grid solutions are plotted above fine-grid solutions for clarity
    for kk = tSteps
        plot(par.xC,C_DM(:,kk),'o','Color',cols(3,:),'MarkerFaceColor',cols(3,:),'LineWidth',10,'MarkerSize',14)
    end

    xlim([0,par.xL])
    ylim([-0.01,1.01])
    set(gca,'FontSize',38,'TickLabelInterpreter','latex','Layer','top')
    xlabel('$x$','Interpreter','latex','FontSize',44)
    ylabel('$c_i^{(k)}$, $\widetilde{c}_i^{(k)}$, $C_m^{(k)}$','Interpreter','latex','FontSize',46)
    set(gcf,'position',[25,25,1600,900],'Renderer','painters')

    if saveFigs
        print(['Case',num2str(testCase),'.eps'],'-depsc2','-r300','-vector')
    end

end

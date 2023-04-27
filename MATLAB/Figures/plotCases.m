% This code will plot the test cases in section 4 [Figure 2] and is reliant
% on the files 'case_i_approx.mat', for i = 1, 2, ..., 6, which contain
% approximate solutions to the 1D ADR problem [Equations 1-4]. If not 
% present, these files can be generated using the code in 
% globalErrorPlot.m.

close all
clearvars
commandwindow; clc;

par.x0 = 0;
par.xL = 30;

par.N = 3601;
par.xF = linspace(par.x0,par.xL,par.N);

par.M = 11;
par.xC = linspace(par.x0,par.xL,par.M);

par.t0 = 0;
par.tK = 18;
par.K = 3600;

cases = [1,2,3,4,5,6];

saveFigs = false;

cols = [
    0.6235    0.6549    0.6784
    0.1373    0.1529    0.1608
    0.0863    0.4431    0.8510
    0.7098    0.0353    0.0471
    ];

for ii = 1:length(cases)

    testCase = cases(ii);

    load(['case_',num2str(testCase),'_approx.mat'])

    run(fullfile(cd,'..','\testCases.m'));

    figure;
    hold on
    box on

    % Fine grid solution plots, from both FVM and dual-grid (DM) methods
    for kk = 1:length(tSteps)
        plot(par.xF,c_FVM(:,tSteps(kk)),'Color',cols(1,:),'LineWidth',10)
        plot(par.xF,c_DM(:,tSteps(kk)),'--','Color',cols(2,:),'LineWidth',10)
    end

    % Coarse-grid solutions are plotted above fine-grid solutions for clarity
    for kk = 1:length(tSteps)
        plot(par.xC,C_DM(:,tSteps(kk)),'o','Color',cols(3,:),'MarkerFaceColor',cols(3,:),'LineWidth',10,'MarkerSize',14)
    end

    xlim([par.x0,par.xL])
    ylim([-0.01,1.01])
    set(gca,'FontSize',38,'TickLabelInterpreter','latex','Layer','top')
    xlabel('$x$','Interpreter','latex','FontSize',44)
    ylabel('$c_i^{(k)}$, $\widetilde{c}_i^{(k)}$, $C_m^{(k)}$','Interpreter','latex','FontSize',46)
    set(gcf,'position',[25,25,1600,900],'Renderer','painters')

    if saveFigs
        print(['Case',num2str(testCase),'.eps'],'-depsc2','-r300','-vector')
    end

end

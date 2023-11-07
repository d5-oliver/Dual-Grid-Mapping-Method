% Note: This file is intended to be executed from generateResults.m.
%
% This code is used to generate maximum absolute difference comparisons
% [Figure 3].

tic

cases = [1,2,3,4,5,6,7];

par.xL = 30;

par.N = 3601;
par.xF = linspace(0,par.xL,par.N);

par.M = 11;
par.xC = linspace(0,par.xL,par.M);

par.tK = 18;
par.K = 3600;
tau = par.tK/par.K;

par.sigma = 1e+10;
par.omega = 1;

tic
dummy = 10/3 + 5/7;
dummy = toc;

if generateNew

    testVals = [3,4,5,6,7,11,13,16,17,19,21,25,26,31,37,41,46,49,51,61,73,76,81,91,101,121,145,151,181,201,226,241,301,361,401,451,601,721,901,1201,1801];

    U = length(testVals);

    V = length(cases);

    errorsDM = zeros(V,U);
    timersDM = zeros(V,U);

    errorsIM = zeros(V,U);
    timersIM = zeros(V,U);

    timersFVM = zeros(V,U);

    for ii = 1:V

        % Generate parameter set for chosen test case
        testCase = cases(ii);
        run([topFolder,'\','testCases.m'])

        tic
        dummy = 10/3 + 5/7;
        dummy = toc;

        tic
        c_FVM = FVM(par);
        timerFVM = toc;

        for jj = 1:U

            % Set current no. of coarse-grid nodes
            par.M = testVals(jj);

            % Define coarse grid (xC)
            par.xC = linspace(0,par.xL,par.M);

            % Set location of coarse-grid nodes in fine grid.
            xSteps = ceil(par.N * par.xC / par.xL);
            xSteps(1) = 1;

            tic
            dummy = 10/3 + 5/7;
            dummy = toc;

            tic
            [C_IM,c_IM] = IM(par);
            timerIM = toc;

            tic
            dummy = 10/3 + 5/7;
            dummy = toc;

            tic
            [C_DM,c_DM] = DM(par);
            timerDM = toc;

            coarseErrorDM = max(abs(c_FVM(xSteps,:) - C_DM),[],2);
            globalErrorDM = max(max(abs(c_FVM - c_DM)));

            coarseErrorIM = max(abs(c_FVM(xSteps,:) - C_IM),[],2);
            globalErrorIM = max(max(abs(c_FVM - c_IM)));

            if par.M == 11

                catErrorDM = [coarseErrorDM; globalErrorDM];
                catErrorIM = [coarseErrorIM;globalErrorIM];

                save([bottomFolder,'\','case_',num2str(testCase),'_approx.mat'],'c_FVM','C_DM','c_DM','C_IM','c_IM')
                save([bottomFolder,'\','case_',num2str(testCase),'_abs_differences.mat'],'catErrorDM','catErrorIM')

            end

            errorsDM(ii,jj) = globalErrorDM;
            timersDM(ii,jj) = timerDM;

            errorsIM(ii,jj) = globalErrorIM;
            timersIM(ii,jj) = timerIM;

            timersFVM(ii,jj) = timerFVM;

        end

        disp(['Test case ',num2str(ii),' complete.'])

    end

    save([bottomFolder,'\','errorsTimers.mat'],'testVals','cases','errorsDM','timersDM','errorsIM','timersIM','timersFVM')

else

    load([bottomFolder,'\','errorsTimers.mat'])

    U = length(testVals);
    V = length(cases);

end

%% Set plot colours

cols = [
    0.6235    0.6549    0.6784
    0.1373    0.1529    0.1608
    0.0863    0.4431    0.8510
    0.7098    0.0353    0.0471
    ];

%% Define error and timer values from local files or freshly generated data

avgErrorValsDM = mean(errorsDM(1:end-1,:));
errorQuantilesDM = quantile(errorsDM(1:end-1,:),[0,1]);
avgTimerValsDM = mean(timersDM(1:end-1,:));
timerQuantilesDM = quantile(timersDM(1:end-1,:),[0,1]);
[minTimer,minIdx] = min(avgTimerValsDM);
minNodes = testVals(minIdx);
minError = avgErrorValsDM(minIdx);

avgErrorValsIM = mean(errorsIM(1:end-1,:));
errorQuantilesIM = quantile(errorsIM(1:end-1,:),[0,1]);
avgTimerValsIM = mean(timersIM(1:end-1,:));
timerQuantilesIM = quantile(timersIM(1:end-1,:),[0,1]);

avgTimerValsFVM = mean(timersFVM(1:end-1,:));
timerQuantilesFVM = quantile(timersFVM(1:end-1,:),[0,1]);

%% Generate maximum abs. difference plot

figure
hold on
box on

fill([testVals,testVals(end:-1:1)],[errorQuantilesIM(1,:),errorQuantilesIM(2,end:-1:1)],cols(4,:),'FaceAlpha',0.25,'EdgeColor','none')
fill([testVals,testVals(end:-1:1)],[errorQuantilesDM(1,:),errorQuantilesDM(2,end:-1:1)],cols(3,:),'FaceAlpha',0.25,'EdgeColor','none')

plot(testVals,avgErrorValsIM,'Color',cols(4,:),'LineWidth',10)
plot(testVals,avgErrorValsDM,'Color',cols(3,:),'LineWidth',10)

plot(minNodes,minError,'o','Color',cols(2,:),'MarkerSize',16,'MarkerFaceColor',cols(2,:))

tickValues = (par.N - 1) ./ ([1801,901,601,361,181,121,61,31,19,9,5,3] - 1) + 1;

xlim([min(testVals),max(testVals)])
ylim([1e-8,1e-0])
set(gca,'YScale','log','YTick',[1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e-0,1e+1],'XScale','log','XTick',tickValues,'XTickLabel',tickValues,'XMinorTick','off','TickLabelInterpreter','latex','FontSize',38,'fontweight','bold','Layer','top','Clipping','on','ClippingStyle','rectangle')
ylabel('$\max\limits_{i,k}|c_i^{(k)} - \widetilde{c}_i^{(k)}|$','Interpreter','latex','FontSize',46,'fontweight','bold')
xlabel('$M$','Interpreter','latex','FontSize',44,'fontweight','bold')

set(gcf,'position',[25,25,1600,900],'Renderer','painters')

xtickangle(0)

if saveFigs

    print([midFolder,'\','Errors.eps'],'-depsc2','-r300','-vector')

end

%% Generate runtime plot

figure
hold on
box on

%fill([testVals,testVals(end:-1:1)],[timerQuantilesFVM(1,:),timerQuantilesFVM(2,end:-1:1)],cols(1,:),'FaceAlpha',0.25,'EdgeColor','none')
fill([testVals,testVals(end:-1:1)],[timerQuantilesIM(1,:),timerQuantilesIM(2,end:-1:1)],cols(4,:),'FaceAlpha',0.25,'EdgeColor','none')
fill([testVals,testVals(end:-1:1)],[timerQuantilesDM(1,:),timerQuantilesDM(2,end:-1:1)],cols(3,:),'FaceAlpha',0.25,'EdgeColor','none')

plot(testVals,avgTimerValsFVM,'Color',cols(1,:),'LineWidth',10)

plot(testVals,avgTimerValsIM,'Color',cols(4,:),'LineWidth',10)

plot(testVals,avgTimerValsDM,'Color',cols(3,:),'LineWidth',10)

plot(minNodes,minTimer,'o','Color',cols(2,:),'MarkerSize',16,'MarkerFaceColor',cols(2,:))

xlim([min(testVals),max(testVals)])
ylim([0,400])
set(gca,'YTick',0:50:400,'XScale','log','XTick',tickValues,'XTickLabel',tickValues,'XMinorTick','off','FontSize',38,'TickLabelInterpreter','latex','fontweight','bold','Layer','top','Clipping','on','ClippingStyle','rectangle')
ylabel('Runtime [s]','Interpreter','latex','FontSize',46,'fontweight','bold')
xlabel('$M$','Interpreter','latex','FontSize',44,'fontweight','bold')

set(gcf,'position',[25,25,1600,900],'Renderer','painters')

xtickangle(0)

if saveFigs

    print([midFolder,'\','Runtimes.eps'],'-depsc2','-r300','-vector')

end

%% Generate legends

figure;
hold on
box on
plot(0,0,'Color',cols(1,:),'linewidth',8)
plot(0,0,'Color',cols(2,:),'linewidth',8)
plot(0,0,'o','Color',cols(3,:),'MarkerFaceColor',cols(3,:),'LineWidth',1.5,'markersize',10)
leg1 = legend('Fine-grid solution$\quad\quad$','Fine-grid DM solution$\quad\quad$','Coarse-grid DM solution','Orientation','horizontal','fontsize',20,'interpreter','latex','fontweight','bold');
set(gca,'XTick',[],'YTick',[],'Layer','top','Clipping','on','ClippingStyle','rectangle')
set(gcf,'Position',get(leg1,'Position').*[0,0,1,1].*get(gcf,'Position'))
set(leg1,'Position',[0,0,1,1])
set(gcf,'Position',get(gcf,'Position')+[100,100,0,0],'Renderer','painters')

if saveFigs

    print([midFolder,'\','Legend_Plots.eps'],'-depsc2','-r300','-vector')

end

figure;
hold on
box off
plot(0,0,'Color',cols(1,:),'linewidth',8)
plot(0,0,'Color',cols(4,:),'linewidth',8)
plot(0,0,'Color',cols(3,:),'LineWidth',8)
leg2 = legend('Fine-grid solutions$\quad\quad$','IM solutions$\quad\quad$','DM solutions','Orientation','horizontal','fontsize',20,'interpreter','latex','fontweight','bold');
set(gca,'XTick',[],'YTick',[],'Layer','top','Clipping','on','ClippingStyle','rectangle')
set(gcf,'Position',get(leg2,'Position').*[0,0,1,1].*get(gcf,'Position'))
set(leg2,'Position',[0,0,1,1])
set(gcf,'Position',get(gcf,'Position')+[100,100,0,0],'Renderer','painters')

if saveFigs

    print([midFolder,'\','Legend_Runtimes.eps'],'-depsc2','-r300','-vector')

end

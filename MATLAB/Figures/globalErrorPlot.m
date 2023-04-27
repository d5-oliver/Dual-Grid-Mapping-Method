% This code is used to generate maximum absolute difference comparisons 
% [Figure 3]. To generate new results, comment lines 4-14 of the runAll.m
% script and set generateNew = true. Otherwise, to simply generate plots,
% set generateNew = false. To also save the generated figures, set 
% saveFigs = true.

clearvars
close all
clc

saveFigs = false;
generateNew = false;

par.N = 3601;

if generateNew

    testVals = [3,4,5,6,7,11,13,16,19,21,25,26,31,37,41,46,49,51,73,76,...
                81,91,101,121,145,151,181,201,226,301,361,401,451,601,...
                721,901,1201,1801];

    U = length(testVals);

    testCaseVals = [1,2,3,4,5,6];
    V = length(testCaseVals);

    errorsDM = zeros(V,U);
    timersDM = zeros(V,U);

    errorsIM = zeros(V,U);
    timersIM = zeros(V,U);

    timersFVM = zeros(V,U);

    for ii = 1:V

        testCase = testCaseVals(ii);

        for jj = 1:U

            par.M = testVals(jj);

            run(fullfile(cd,'..','\runAll.m'))

            xSteps = ceil(par.N * par.xC / (par.xL - par.x0));
            xSteps(1) = 1;

            coarseErrorDM = max(abs(c_FVM(xSteps,tSteps) - C_DM(:,tSteps)),[],2);
            globalErrorDM = max(max(abs(c_FVM(:,tSteps) - c_DM(:,tSteps))));

            coarseErrorIM = max(abs(c_FVM(xSteps,tSteps) - C_IM(:,tSteps)),[],2);
            globalErrorIM = max(max(abs(c_FVM(:,tSteps) - c_IM(:,tSteps))));

            if par.M == 11

                catErrorDM = [coarseErrorDM; globalErrorDM];
                catErrorIM = [coarseErrorIM;globalErrorIM];

                save(['case_',num2str(testCase),'_approx.mat'],'c_FVM','C_DM','c_DM','C_IM','c_IM')
                save(['case_',num2str(testCase),'_abs_differences.mat'],'catErrorDM','catErrorIM')

            end

            errorsDM(ii,jj) = globalErrorDM;
            timersDM(ii,jj) = timerLocal;

            errorsIM(ii,jj) = globalErrorIM;
            timersIM(ii,jj) = timerInterp;

            timersFVM(ii,jj) = timerFVM;

        end

        disp(['Test case ',num2str(ii),' complete.'])

    end

    save('errorsTimers.mat','testVals','testCaseVals','errorsDM','timersDM','errorsIM','timersIM','timersFVM')

else

    load('errorsTimers.mat')

    U = length(testVals);
    V = length(testCaseVals);

end

%% Set plot colours

cols = [
        0.6235    0.6549    0.6784
        0.1373    0.1529    0.1608
        0.0863    0.4431    0.8510
        0.7098    0.0353    0.0471
       ];

%% Define error and timer values from local files or freshly generated data

avgErrorValsDM = mean(errorsDM);
errorQuantilesDM = quantile(errorsDM,[0,1]);
avgTimerValsDM = mean(timersDM);
timerQuantilesDM = quantile(timersDM,[0,1]);
[minTimer,minIdx] = min(avgTimerValsDM);
minNodes = testVals(minIdx);
minError = avgErrorValsDM(minIdx);

avgErrorValsIM = mean(errorsIM);
errorQuantilesIM = quantile(errorsIM,[0,1]);
avgTimerValsIM = mean(timersIM);
timerQuantilesIM = quantile(timersIM,[0,1]);

avgTimerValsFVM = mean(timersFVM);
timerQuantilesFVM = quantile(timersFVM,[0,1]);

%% Generate maximum abs. difference plot

figure
hold on
box on

fill([testVals,testVals(end:-1:1)],[errorQuantilesIM(1,:),errorQuantilesIM(2,end:-1:1)],cols(4,:),'FaceAlpha',0.25,'EdgeColor','none')
fill([testVals,testVals(end:-1:1)],[errorQuantilesDM(1,:),errorQuantilesDM(2,end:-1:1)],cols(3,:),'FaceAlpha',0.25,'EdgeColor','none')

plot(testVals,avgErrorValsIM,'Color',cols(4,:),'LineWidth',10)
plot(testVals,avgErrorValsDM,'Color',cols(3,:),'LineWidth',10)

plot(minNodes,minError,'o','Color',cols(2,:),'MarkerSize',16,'MarkerFaceColor',cols(2,:))

tickValues = (par.N-1) ./ [1800,900,600,360,180,120,60,30,18,8,4,2] + 1;

xlim([min(testVals),max(testVals)])
ylim([1e-8,1e-0])
set(gca,'YScale','log','YTick',[1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e-0,1e+1],'XScale','log','XTick',tickValues,'XTickLabel',tickValues,'XMinorTick','off','TickLabelInterpreter','latex','FontSize',38,'fontweight','bold','Layer','top','Clipping','on','ClippingStyle','rectangle')
ylabel('$\max\limits_{i,k}|c_i^{(k)} - \widetilde{c}_i^{(k)}|$','Interpreter','latex','FontSize',46,'fontweight','bold')
xlabel('$M$','Interpreter','latex','FontSize',44,'fontweight','bold')

set(gcf,'position',[25,25,1600,900],'Renderer','painters')

xtickangle(0)

if saveFigs

    print('Errors.eps','-depsc2','-r300','-vector')

end

%% Generate runtime plot

figure
hold on
box on

fill([testVals,testVals(end:-1:1)],[timerQuantilesFVM(1,:),timerQuantilesFVM(2,end:-1:1)],cols(1,:),'FaceAlpha',0.25,'EdgeColor','none')
fill([testVals,testVals(end:-1:1)],[timerQuantilesIM(1,:),timerQuantilesIM(2,end:-1:1)],cols(4,:),'FaceAlpha',0.25,'EdgeColor','none')
fill([testVals,testVals(end:-1:1)],[timerQuantilesDM(1,:),timerQuantilesDM(2,end:-1:1)],cols(3,:),'FaceAlpha',0.25,'EdgeColor','none')

plot(testVals,avgTimerValsFVM,'Color',cols(1,:),'LineWidth',10)

plot(testVals,avgTimerValsIM,'Color',cols(4,:),'LineWidth',10)

plot(testVals,avgTimerValsDM,'Color',cols(3,:),'LineWidth',10)

plot(minNodes,minTimer,'o','Color',cols(2,:),'MarkerSize',16,'MarkerFaceColor',cols(2,:))

xlim([min(testVals),max(testVals)])
ylim([0,325])
set(gca,'YTick',0:25:325,'XScale','log','XTick',tickValues,'XTickLabel',tickValues,'XMinorTick','off','FontSize',38,'TickLabelInterpreter','latex','fontweight','bold','Layer','top','Clipping','on','ClippingStyle','rectangle')
ylabel('Runtime [s]','Interpreter','latex','FontSize',46,'fontweight','bold')
xlabel('$M$','Interpreter','latex','FontSize',44,'fontweight','bold')

set(gcf,'position',[25,25,1600,900],'Renderer','painters')

xtickangle(0)

if saveFigs

    print('Runtimes.eps','-depsc2','-r300','-vector')

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

    print('Legend_Plots.eps','-depsc2','-r300','-vector')

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

    print('Legend_Runtimes.eps','-depsc2','-r300','-vector')

end

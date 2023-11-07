% Note: This file is intended to be executed from generateResults.m.
%
% Generates a latex table (requiring the tabularx package) of solution
% values at select time steps and at each of the coarse-grid nodes. The
% table assumes the coarse grid is defined with 11 nodes over the same
% range as the test cases in the paper, 0 < x < 30. This code is reliant on
% the files 'case_j_approx.mat', for j = 1, 2, ..., 6, which contains
% approximation values obtained using the dual-grid mapping method (IM and
% DM) and the fine-grid solution method. If not present, these files can be
% generated using the code in globalErrorPlot.m.

cases = [1,2,3,4,5,6];%,7];

par.xL = 30;

par.N = 3601;
par.xF = linspace(0,par.xL,par.N);

par.M = 11;
par.xC = linspace(0,par.xL,par.M);

par.tK = 18;
par.K = 3600;

par.sigma = 1e+10;
par.omega = 1;

xSteps = find(mod(par.xF,1) == 0);
xCol = par.xF(xSteps);

% Finding the harmonic average of the diffusivities for each test problem
aveD = zeros(1,length(cases));
L = par.xL;
h = par.xF(2) - par.xF(1);
for testCase = cases
    run([topFolder,'\','testCases.m'])
    D = par.D(par.xF);
    aveD(testCase) = L / (h/2 * (1./D(1) + 2*sum(1./D(2:end-1)) + 1./D(end)));
end
% Note: we round the time value to 4 decimal places so that we can use the
% same (reasonably sized) time step in each simulation.
% Using Fourier number formula to convert dimensionless time to a standard
% time value.
tVals = round(L^2 * 0.005 ./ aveD,4);
% Setting number of time steps so that tau=0.0001 in each test problem
% (round(*) simply converts to integer in this case).
KVals = round(tVals / 0.0001);

U = length(xSteps) + 2;
V = 2*length(cases) + 1;

tableCoarse = cell(U,V);

tableCoarse{2,1} = '$x$';

for jj = 1:length(cases)

    testCase = cases(jj);
    run([topFolder,'\','testCases.m'])

    par.tK = tVals(jj);
    par.K = KVals(jj);

    if generateNew
        c_FVM = FVM_Sparse(par);
        [C_DM,c_DM] = DM_Sparse(par);

        save(['case_',num2str(cases(jj)),'_dimensionless.mat'],'c_FVM','C_DM','c_DM')
    else
        load(['case_',num2str(cases(jj)),'_dimensionless.mat'])
    end

    tableCoarse{1,2*(jj-1)+2} = ['\\multicolumn{2}{l}{Case ',num2str(testCase),'}'];
    tableCoarse{2,2*(jj-1)+2} = 'FGS';
    tableCoarse{2,2*(jj-1)+3} = 'DM';

    for ii = 1:length(xCol)

        if mod(xCol(ii),3) == 0
            tableCoarse{ii+2,1} = ['\\bf{',num2str(xCol(ii)),'}'];
        else
            tableCoarse{ii+2,1} = num2str(xCol(ii));
        end

        tmp = sprintf('%.4f',c_FVM(xSteps(ii),end));
        tableCoarse{ii+2,2*(jj-1)+2} = ['$',tmp,'$'];

        tmp = sprintf('%.4f',c_DM(xSteps(ii),end));
        tableCoarse{ii+2,2*(jj-1)+3} = ['$',tmp,'$'];

    end

end

latexCoarse = fopen([midFolder,'\','tableCoarse.tex'],'w');

fprintf(latexCoarse,['\\begin{tabularx}{\\linewidth}{',repmat('X',1,V),'}\n\\toprule\n & ']);

for jj = 1:length(cases)
    if jj == length(cases)
        fprintf(latexCoarse,[tableCoarse{1,2*(jj-1)+2},' ']);
    else
        fprintf(latexCoarse,[tableCoarse{1,2*(jj-1)+2},' & ']);
    end
end
fprintf(latexCoarse,'\\\\\n');

for jj = 1:length(cases)
    fprintf(latexCoarse,['\\cmidrule(lr){',num2str(2*(jj-1)+2),'-',num2str(2*(jj-1)+3),'} ']);
end
fprintf(latexCoarse,'\n');

for jj = 1:V
    if jj == V
        fprintf(latexCoarse,[tableCoarse{2,jj},' ']);
    else
        fprintf(latexCoarse,[tableCoarse{2,jj},' & ']);
    end
end
fprintf(latexCoarse,'\\\\\n\\midrule\n');

for ii = 3:U
    for jj = 1:V
        if jj == V
            fprintf(latexCoarse,[tableCoarse{ii,jj},' \\\\\n']);
        else
            fprintf(latexCoarse,[tableCoarse{ii,jj},' & ']);
        end
    end
end
fprintf(latexCoarse,'\\bottomrule\n\\end{tabularx}\n\n');

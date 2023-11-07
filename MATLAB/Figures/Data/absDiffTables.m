% Note: This file is intended to be executed from generateResults.m.
%
% Generates a latex table (requiring the tabularx package) of maximum
% absolute differences over all time steps at each of the coarse-grid nodes
% (and also the overall maximum absolute difference between the benchmark
% and reconstructed fine-scale solution for each separate case, which is
% stored in the 'Global' column). The difference table assumes the coarse
% grid is defined with 11 nodes over the same range as the test cases in
% the paper, 0 < x < 30. This code is reliant on the files
% 'case_j_abs_differences.mat', for j = 1, 2, ..., 6, which contains
% maximum absolute differences between approximations obtained using the
% dual-grid mapping method (IM and DM) and the fine-grid solution method.
% If not present, these files can be generated using the code in
% globalErrorPlot.m.

% Set test problem parameters
cases = [1,2,3,4,5,6,7];

par.xL = 30;

par.N = 3601;
par.xF = linspace(0,par.xL,par.N);

par.M = 11;
par.xC = linspace(0,par.xL,par.M);

par.tK = 18;
par.K = 3600;

par.sigma = 1e+10;
par.omega = 1;

xCol = par.xC;

testCase = 1;
run([topFolder,'\','testCases.m'])

U = length(xCol) + 2;
V = length(cases) + 1;

diffTableDM = cell(U,V);
diffTableIM = cell(U,V);

diffTableDM{1,1} = 'Case';
diffTableIM{1,1} = 'Case';

for jj = 2:V
    tmp = num2str(cases(jj-1));
    diffTableDM{1,jj} = tmp;
    diffTableIM{1,jj} = tmp;
end

for ii = 2:U-1
    tmp = num2str(xCol(ii-1));
    diffTableDM{ii,1} = tmp;
    diffTableIM{ii,1} = tmp;
end

diffTableDM{end,1} = 'Global';
diffTableIM{end,1} = 'Global';

for jj = 2:V

    load([bottomFolder,'\','case_',num2str(jj-1),'_abs_differences.mat'])

    for ii = 2:U

        strValDM = sprintf('%.2d',catErrorDM(ii-1));
        diffTableDM{ii,jj} = [strValDM(1:4),' \\times 10^{',num2str(str2double(strValDM(end-2:end))),'}'];

        strValIM = sprintf('%.2d',catErrorIM(ii-1));
        diffTableIM{ii,jj} = [strValIM(1:4),' \\times 10^{',num2str(str2double(strValIM(end-2:end))),'}'];

    end

end

latexStringDM = fopen([midFolder,'\','absDiffTable_DM.tex'],'w');
fprintf(latexStringDM,['\\begin{tabularx}{\\linewidth}{',repmat('X',1,V+1),'} \n \\toprule \n ']);

latexStringIM = fopen([midFolder,'\','absDiffTable_IM.tex'],'w');
fprintf(latexStringIM,['\\begin{tabularx}{\\linewidth}{',repmat('X',1,V+1),'} \n \\toprule \n ']);

for jj = 1:V
    if (jj > 1) && (jj < V)
        fprintf(latexStringDM,[num2str(diffTableDM{1,jj}),' & ']);
        fprintf(latexStringIM,[num2str(diffTableIM{1,jj}),' & ']);
    elseif (jj == V)
        fprintf(latexStringDM,num2str(diffTableDM{1,jj}));
        fprintf(latexStringIM,num2str(diffTableIM{1,jj}));
    else
        fprintf(latexStringDM,[diffTableDM{1,jj},' & ']);
        fprintf(latexStringIM,[diffTableIM{1,jj},' & ']);
    end
end

fprintf(latexStringDM,' \\\\ \n \\midrule \n ');
fprintf(latexStringIM,' \\\\ \n \\midrule \n ');

for ii = 2:U
    for jj = 1:V
        if (jj == 1) && (ii ~= U)
            fprintf(latexStringDM,['$x = ',diffTableDM{ii,jj},'$ & ']);
            fprintf(latexStringIM,['$x = ',diffTableIM{ii,jj},'$ & ']);
        elseif (jj == 1) && (ii == U)
            fprintf(latexStringDM,[diffTableDM{ii,jj},' & ']);
            fprintf(latexStringIM,[diffTableIM{ii,jj},' & ']);
        elseif (jj == V)
            fprintf(latexStringDM,['$',diffTableDM{ii,jj},'$ ']);
            fprintf(latexStringIM,['$',diffTableIM{ii,jj},'$ ']);
        else
            fprintf(latexStringDM,['$',diffTableDM{ii,jj},'$ & ']);
            fprintf(latexStringIM,['$',diffTableIM{ii,jj},'$ & ']);
        end
    end
    fprintf(latexStringDM,'\\\\ \n ');
    fprintf(latexStringIM,'\\\\ \n ');
end

fprintf(latexStringDM,'\\bottomrule \n ');
fprintf(latexStringIM,'\\bottomrule \n ');

fprintf(latexStringDM,'\\end{tabularx}');
fprintf(latexStringIM,'\\end{tabularx}');

fclose(latexStringDM);
fclose(latexStringIM);

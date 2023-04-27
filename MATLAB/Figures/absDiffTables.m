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

U = 12;
V = 6;

diffTableDM = cell(U,V);
diffTableIM = cell(U,V);
diff0 = linspace(0,30,11);

for jj = 1:length(diff0)

    diffTableDM{jj,1} = num2str(diff0(jj));
    diffTableIM{jj,1} = num2str(diff0(jj));

end

diffTableDM{end,1} = 'Global';
diffTableIM{end,1} = 'Global';

for jj = 1:V

    load(['case_',num2str(jj),'_abs_differences.mat'])

    for ii = 1:U

        strValDM = sprintf('%.2d',catErrorDM(ii));
        diffTableDM{ii,jj+1} = [strValDM(1:4),' \times 10^{-',num2str(str2double(strValDM(end-1:end))),'}'];

        strValIM = sprintf('%.2d',catErrorIM(ii));
        diffTableIM{ii,jj+1} = [strValIM(1:4),' \times 10^{-',num2str(str2double(strValIM(end-1:end))),'}'];

    end

end

diffTableDM = ['Case',num2cell(1:6);diffTableDM];
latexStringDM = ['\begin{tabularx}{\linewidth}{',repmat('X',1,V+1),'} \toprule '];

diffTableIM = ['Case',num2cell(1:6);diffTableIM];
latexStringIM = ['\begin{tabularx}{\linewidth}{',repmat('X',1,V+1),'} \toprule '];

for ii = 1:U+1

    for jj = 1:V+1

        if isnumeric(diffTableDM{ii,jj})

            diffTableDM{ii,jj} = num2str(diffTableDM{ii,jj});

            diffTableIM{ii,jj} = num2str(diffTableIM{ii,jj});

        end

        if (ii == 1) && (jj == 1)

            latexStringDM = append(latexStringDM,[diffTableDM{ii,jj},' & ']);

            latexStringIM = append(latexStringIM,[diffTableIM{ii,jj},' & ']);

        elseif (ii == 1) && (jj < V+1)

            latexStringDM = append(latexStringDM,['$',diffTableDM{ii,jj},'$ & ']);

            latexStringIM = append(latexStringIM,['$',diffTableIM{ii,jj},'$ & ']);

        elseif (ii == 1) && (jj == V+1)

            latexStringDM = append(latexStringDM,['$',diffTableDM{ii,jj},'$ \\ \midrule ']);

            latexStringIM = append(latexStringIM,['$',diffTableIM{ii,jj},'$ \\ \midrule ']);

        elseif (ii == U+1) && (jj == 1)

            latexStringDM = append(latexStringDM,[diffTableDM{ii,jj},' & ']);

            latexStringIM = append(latexStringIM,[diffTableIM{ii,jj},' & ']);

        elseif (ii == U+1) && (jj == V+1)

            latexStringDM = append(latexStringDM,['$',diffTableDM{ii,jj},'$ \\ \bottomrule']);

            latexStringIM = append(latexStringIM,['$',diffTableIM{ii,jj},'$ \\ \bottomrule']);

        elseif jj == V+1

            latexStringDM = append(latexStringDM,['$',diffTableDM{ii,jj},'$ \\ ']);

            latexStringIM = append(latexStringIM,['$',diffTableIM{ii,jj},'$ \\ ']);

        elseif jj == 1

            latexStringDM = append(latexStringDM,['$x = ',diffTableDM{ii,jj},'$ & ']);

            latexStringIM = append(latexStringIM,['$x = ',diffTableIM{ii,jj},'$ & ']);

        else

            latexStringDM = append(latexStringDM,['$',diffTableDM{ii,jj},'$ & ']);

            latexStringIM = append(latexStringIM,['$',diffTableIM{ii,jj},'$ & ']);

        end

    end

end

latexStringDM = append(latexStringDM,' \end{tabularx}');
latexStringIM = append(latexStringIM,' \end{tabularx}');

disp(latexStringDM)
disp(latexStringIM)

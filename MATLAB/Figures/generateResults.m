% Generate figures and tables used within the paper. See individual files
% for more information. Comment lines 1-21 of runAll.m in the parent folder
% before running.

clearvars
close all
clc

% Set folder locations based on current directory
midFolder = cd;
topFolder = midFolder(1:end-8);
bottomFolder = [midFolder,'\Data'];

% Choose to generate all new data or generate results from pre-existing
% data files. Note: if there are no .mat files present in the Data folder,
% new data will need to be generated (WARNING: time and resource
% intensive).
generateNew = 0;

% Choose to save figures in .eps format.
saveFigs = 1;

disp('Starting max. abs. difference test...')
run([bottomFolder,'\','globalErrorPlot.m'])
disp('Max. abs. difference test complete.')

disp('Generating test case plots...')
run([bottomFolder,'\','plotCases.m'])
disp('Test case plotting complete.')

disp('Starting order of accuracy test...')
run([bottomFolder,'\','orderTest.m'])
disp('Order of accuracy test complete.')

disp('Generating max. abs. difference tables...')
run([bottomFolder,'\','absDiffTables.m'])
disp('Max. abs. difference table generation complete.')

disp('Generating test case solution value tables...')
run([bottomFolder,'\','approxTables.m'])
disp('Solution value table generation complete.')

disp('All tests complete.')

function R = project_AddPaths(R)
%
% first obtain name of the current user
U = char(java.lang.System.getProperty('user.name'));
switch U
    case 'DESKTOP-1QJTIMO'
        gitpath =  'C:\Users\Tim West\Documents\GitHub';
        madpath = 'C:\Users\Tim West\Documents\MATLAB ADDONS';
        spmpath = 'C:\Users\Tim West\Documents\GitHub\spm12';
    case 'DESKTOP-0HO6J14'
        gitpath =  'D:\GITHUB';
        madpath = 'C:\Users\timot\OneDrive\Documents\Work\MATLAB ADDONS';
        spmpath = 'D:\GITHUB\spm12';
    case 'jyao'                                                             % for all of jyao's scripts
        gitpath = fullfile('~', 'local', 'gitprojects');
        spmpath = fullfile(gitpath, 'spm12');
        
    otherwise
        error('You need to add your PC specific paths - ...please look at project_AddPaths.m for template')
end
% Sometimes the main project is some other folder to Github default
if ~exist('gitpath2','var')
    gitpath2 = gitpath;
end

R.rootn = fullfile(gitpath2, 'StimulatingAtTheRightTime_West2022');

% Add the root
addpath(genpath(R.rootn))

% Add SPM/Fieldtrip
pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(spmpath, pathCell));
if ~onPath; addpath(spmpath); spm eeg; close all; end


function [R, m, p] = loadModelForDRL(Rorg)
%% Load Model Parameters
fittedModelPath = fullfile(Rorg.rootn, 'data', 'modelfit', ...
    'SimModelData_M10.mat');
load(fittedModelPath,'R','m','p')
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;

% Simulate the Base Model
R.obs.trans.norm = 1; % Normalize the output spectra
R.obs.trans.gauss = 1; % Smooth
R.obs.obsstates = [1:6]; % All 6 nodes are observed
R.chloc_name = R.chsim_name; % Ensure sim names match to output names
R.obs.gainmeth = {'unitvar'};

end
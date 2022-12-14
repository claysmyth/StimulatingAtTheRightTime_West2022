%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

fprintf('\n>>Initializing models\n')
R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

NUM_REP = 2;
%% Main loop

% simulate with cDBS parameters
% generate around 800 cDBS samples
% each loop gives around 250

% stim amp to play with and make copies of initial structure
vecR = deepCopyStruct(R, NUM_REP);

tic;
parfor i = 1:NUM_REP
    % unpack a little bit
    R_i = vecR{i};
    
    % simulate
    [X_i, R_i] = DRL_sim_genInitEvalBuffer(R_i);

    % pack to outer list
    vecX_eval{i} = X_i;
    vecR{i} = R_i;
end

% combine into single structure
X_eval = combineXStruct(vecX_eval, vecR{1});

% get output path
outputPath = fullfile(R.DRL.outPath, 'initBufferEval');
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

% % save all structs as output
R_eval = vecR{1};
save(fullfile(outputPath, 'X_eval_init'), 'X_eval', '-v7.3');
save(fullfile(outputPath, 'R_eval_init'), 'R_eval', '-v7.3');

% make a separate copy for python
for condsel = 1:numel(R.condnames)
    dataPython{condsel}.prev_xsims_gl.S = cat(3, X_eval{condsel}.xsims_gl.S{:});
    dataPython{condsel}.metadata = X_eval{condsel}.metadata;
end

save(fullfile(outputPath, 'dataPython_eval_init'), 'dataPython', '-v6');
timeElapsed = toc;


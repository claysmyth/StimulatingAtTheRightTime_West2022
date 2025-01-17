%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

fprintf('\n>>Initializing models\n')
R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

%% Main loop

% simulate with cDBS parameters
% generate around 800 cDBS samples
% each loop gives around 250

% stim amp to play with and make copies of initial structure
vecAmpTest = linspace(1e3, 2e3, 3);
vecR = deepCopyStruct(R, numel(vecAmpTest));

tic;
parfor i = 1:numel(vecAmpTest)
    % unpack a little bit
    R_i = vecR{i};
    stimAmp_i = vecAmpTest(i);
    
    % simulate
    [X_i, X_stim_i, A_i, R_i] = DRL_sim_genInitBuffer(R_i, ...
        'SScomb', 1, 'stimAmp', stimAmp_i);

    % pack to outer list
    vecX_cDBS{i} = X_i;
    vecX_stim_cDBS{i} = X_stim_i;
    vecA_cDBS{i} = A_i;
    vecR{i} = R_i;
end

% combine into single structure
X_comb_cDBS = combineXStruct(vecX_cDBS, vecR{1});
X_stim_comb_cDBS = combineXStimStruct(vecX_stim_cDBS, vecR{1});
A_comb_cDBS = combineAStruct(vecA_cDBS, vecR{1});

% simulate with random stim parameters
% generate around 800 samples
vecR = deepCopyStruct(R, numel(vecAmpTest));
parfor i = 1:numel(vecAmpTest)
    % unpack a little bit
    R_i = vecR{i};
    stimAmp_i = vecAmpTest(i);
    
    % simulate
    [X_i, X_stim_i, A_i, R_i] = DRL_sim_genInitBuffer(R_i, ...
        'SScomb', 2, 'stimAmp', stimAmp_i);

    % pack to outer list
    vecX_rDBS{i} = X_i;
    vecX_stim_rDBS{i} = X_stim_i;
    vecA_rDBS{i} = A_i;
    vecR{i} = R_i;
end

% combine into single structure
X_comb_rDBS = combineXStruct(vecX_rDBS, vecR{1});
X_stim_comb_rDBS = combineXStimStruct(vecX_stim_rDBS, vecR{1});
A_comb_rDBS = combineAStruct(vecA_rDBS, vecR{1});

% combine everything into a single last structure
X = combineXStruct({X_comb_cDBS, X_comb_rDBS}, vecR{1});
nextX = combineXStimStruct({X_stim_comb_cDBS, X_stim_comb_rDBS}, vecR{1});
A = combineAStruct({A_comb_cDBS, A_comb_rDBS}, vecR{1});

% get output path
outputPath = fullfile(R.DRL.outPath, 'initBuffer');
if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

% save all structs as output
R = vecR{1};
save(fullfile(outputPath, 'X_init'), 'X', '-v7.3');
save(fullfile(outputPath, 'nextX_init'), 'nextX', '-v7.3');
save(fullfile(outputPath, 'A_init'), 'A', '-v7.3');
save(fullfile(outputPath, 'R_init'), 'R', '-v7.3');

% now compute reward
C = DRL_calc_reward(X, nextX, vecR{1}, 'rewardMethod', 3);

% make a separate copy for python
for condsel = 1:numel(R.condnames)
    dataPython{condsel} = deepCopyStruct(nextX{condsel}, 1);

    dataPython{condsel}.prev_xsims_gl.S = cat(3, X{condsel}.xsims_gl.S{:});
    dataPython{condsel}.metadata = X{condsel}.metadata;

    % also convert cell into array for python
    dataPython{condsel}.xsims_gl.S = cat(3, dataPython{condsel}.xsims_gl.S{:});

    % obtain actions
    dataPython{condsel}.uexs.S = cat(3, A{condsel}.uexs.S{:});

    % obtain rewards
    dataPython{condsel}.reward = C{condsel};
end

save(fullfile(outputPath, 'dataPython_init'), 'dataPython', '-v6');
timeElapsed = toc;


%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

fprintf('\n>>Initializing models\n')
R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

% change this below
rewardMethod = 4;

% parse the reward method
if rewardMethod == 1
    strRewardMethod = 'power';

elseif rewardMethod == 2
    strRewardMethod = 'powerNorm';

elseif rewardMethod == 3
    strRewardMethod = 'burst';

elseif rewardMethod == 4
    strRewardMethod = 'burstNorm';

else
    error('Unknown reward method')
end

%% Load file and recalculate reward

outputPath = fullfile(R.DRL.outPath, 'initBuffer');

A = load(fullfile(outputPath, 'A_init')).A;
R_temp = load(fullfile(outputPath, 'R_init')).R;
X = load(fullfile(outputPath, 'X_init')).X;
nextX = load(fullfile(outputPath, 'nextX_init')).nextX;

% now compute reward
C = DRL_calc_reward(X, nextX, R_temp, 'rewardMethod', rewardMethod);

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

fOutput = sprintf('dataPython_%s_init', strRewardMethod);

save(fullfile(outputPath, fOutput), 'dataPython', '-v6');
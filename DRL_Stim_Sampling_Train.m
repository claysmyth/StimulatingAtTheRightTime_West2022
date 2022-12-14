%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

fprintf('\n>>Initializing models\n')
R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

inputPath = fullfile(R.DRL.outPath, 'model');
initPath = fullfile(R.DRL.outPath, 'initBuffer');

% change these flags
strExpName = 'debug';
rewardMethod = 1;


%% sanity check


R_temp = load(fullfile(initPath, 'R_init.mat')).R;
stimCh = R_temp.IntP.phaseStim.sensStm(2);

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

%% Load all stim patterns

% try to find right directory
strVecTrain = glob(fullfile(inputPath, strExpName, '**', 'A*'));

idxEpoch = [];
upScaleFactor = fix(fix(1/R.IntP.dt) / R.DRL.dsFs);
for i = 1:numel(strVecTrain)
    currStrTrain = strVecTrain{i};
    uexsTemp = load(currStrTrain).uexs;
    
    % first get the index right
    idxEpoch = [idxEpoch, i * ones(1, size(uexsTemp, 1))];

    for condsel = 1:numel(R.condnames)    
        % next make stim pattern
        S = {};
        for j = 1:size(uexsTemp, 1)
            % current stim from model
            uexsTempCurr = uexsTemp(j, :);
            
            % get the indices for interpolation
            idxCurr = 1:upScaleFactor:(upScaleFactor * size(uexsTemp, 2));
            assert(size(idxCurr, 2) == size(uexsTemp, 2));
            idxQuery = 1:(upScaleFactor * size(uexsTemp, 2));
    
            % perform interpolation
            uexsCurrUs = pchip(idxCurr, uexsTempCurr, idxQuery);
    
            % form curren stim
            uexsFull = zeros(size(uexsCurrUs, 2), numel(R.chsim_name));
            uexsFull(:, stimCh) = uexsCurrUs;
    
            S{j} = uexsFull;
        end
        
        % after all indices have been parsed
        A_i{condsel}.uexs.S = S;
    end

    vecA{i} = A_i;
end

% combine loaded A structure
A_comb = combineAStruct(vecA, R_temp);

%% resimulate with new stim

[~, m, p] = loadModelForDRL(R);                                            
X = load(fullfile(initPath, 'X_init.mat')).X;

% resimulate
[X_stim, ~, R] = DRL_resim_bufferStim(X, R_temp, m, p, ...
    'SScomb', 3, 'exStimStruct', A_comb);

% downsample
X_stim = DRL_downSampleData_Redact(X_stim, R);

% recalculate reward
R_new = DRL_calc_reward(X, X_stim, R, 'rewardMethod', rewardMethod, ...
    'boolSave', false, 'init', 0);

%% form summary statistics

fOutput = sprintf('%s_summary.mat', strExpName);
pfeOutput = fullfile(inputPath, strExpName, fOutput);

% if already created then load
if exist(pfeOutput, 'file')
    summary = load(pfeOutput).summary;
end

% compute training loss
R_new_cond = R_new{1};
trainAvgReturn = [];
for i = 1:numel(strVecTrain)
    idxEpochCurr = idxEpoch == i;

    trainAvgReturn = [trainAvgReturn, mean(R_new_cond(idxEpochCurr))];

end

% save summary output
summary.trainAvgReturn = trainAvgReturn;
save(pfeOutput, 'summary');

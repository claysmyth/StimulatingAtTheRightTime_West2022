%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

fprintf('\n>>Initializing models\n')
R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

inputPath = fullfile(R.DRL.outPath, 'model');
initPath = fullfile(R.DRL.outPath, 'initBuffer');
initEvalPath = fullfile(R.DRL.outPath, 'initBufferEval');

% change these flags
useRNN = false;
rewardMethod = 2;

%% sanity check

R_temp = load(fullfile(initPath, 'R_init.mat')).R;
stimCh = R_temp.IntP.phaseStim.sensStm(2);

% get model name
if useRNN
    strModel = 'rnn';
else
    strModel = 'mlp';
end

% get reward name
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


% form correct experiment name
strExpName = sprintf('%s_%s', strModel, strRewardMethod);

%% Load all stim patterns

% try to find right directory
strVecTrain = glob(fullfile(inputPath, strExpName, '**', 'evalA*'));

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
        evalA_i{condsel}.uexs.S = S;
    end

    vec_EvalA{i} = evalA_i;
end


%% resimulate with new stim
parfor i = 1:numel(vec_EvalA)
    % load data again
    [~, m_i, p_i] = loadModelForDRL(R);                                            
    X_eval_i = load(fullfile(initEvalPath, 'X_eval_init.mat')).X_eval;

    % resimulate
    [X_stim_i, ~, R_i] = DRL_resim_bufferStim(X_eval_i, R_temp, m_i, p_i, ...
        'SScomb', 3, 'exStimStruct', vec_EvalA{i});

    % downsample
    X_stim_i = DRL_downSampleData_Redact(X_stim_i, R_i);

    % recalculate reward
    R_new = DRL_calc_reward(X_eval_i, X_stim_i, R_i, 'rewardMethod', rewardMethod, ...
        'boolSave', false, 'init', 0);
    R_new_cond = R_new{1};

    evalAvgReturn{i} = mean(R_new_cond);
end

%% form summary statistics

fOutput = sprintf('%s_summary.mat', strExpName);
pfeOutput = fullfile(inputPath, strExpName, fOutput);

% if already created then load
if exist(pfeOutput, 'file')
    summary = load(pfeOutput).summary;
end

% save summary output
evalAvgReturn = cat(2, evalAvgReturn{:});
summary.evalAvgReturn = evalAvgReturn;
save(pfeOutput, 'summary');

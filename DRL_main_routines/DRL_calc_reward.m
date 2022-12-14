function C = DRL_calc_reward(X, X_stim, R, varargin)
%% Input parsing
% Handle the optional inputs
p = inputParser;
p.KeepUnmatched = true;

addParameter(p, 'rewardMethod', 1, ...                                            % stimulation type
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'boolSave', true, ...                                            % stimulation type
    @(x) validateattributes(x, {'logical'}, {'nonempty'}));

addParameter(p, 'init', 1, ...                                            % stimulation type
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'strSubFolder', '', ...                                            % stimulation type
    @(x) validateattributes(x, {'char'}, {'nonempty'}));

parse(p,varargin{:});

% Handles incorrect inputs
UnmatchedParam = fieldnames(p.Unmatched);
if ~isempty(UnmatchedParam)
    error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
end

% unpacking variable
rewardMethod = p.Results.rewardMethod;
boolSave = p.Results.boolSave;
init = p.Results.init;
strSubFolder = p.Results.strSubFolder;

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


%% Code Body

% loop through the conditions first
for condsel = 1:numel(R.condnames)
    % unpack the struct first
    X_cond_xsims_gl_S = X{condsel}.xsims_gl.S;
    X_cond_xsims_gl_burstThresh = X{condsel}.xsims_gl.burstThresh;
    X_stim_cond_xsims_gl_S = X_stim{condsel}.xsims_gl.S;

    % sanity check
    assert(numel(X_cond_xsims_gl_S) == numel(X_stim_cond_xsims_gl_S));
    numEpoch = numel(X_cond_xsims_gl_S);
    senseSite = R.IntP.phaseStim.sensStm(1);
    fs = R.DRL.dsFs;
    
    % make the filters
    bpFilt_Beta1 = designfilt('bandpassfir','FilterOrder', fix(fs / 10), ...
        'CutoffFrequency1', 13,'CutoffFrequency2', 20, ...
         'SampleRate', fs);
    bpFilt_Beta2 = designfilt('bandpassfir','FilterOrder', fix(fs / 10), ...
        'CutoffFrequency1', 20,'CutoffFrequency2', 30, ...
         'SampleRate', fs);
    
    reward_cond = {};
    for i = 1:numEpoch
        % obtain the relevant data
        currS = X_cond_xsims_gl_S{i}(senseSite, :);
        currS_stim = X_stim_cond_xsims_gl_S{i}(senseSite, :);
        currBurstThresh = X_cond_xsims_gl_burstThresh{i}(senseSite);
        
        % obtain valid index unaffected by filtering
        idxValid = ones(size(currS));
        idxPre = 1:fix(fs / 10);
        idxPost = (numel(idxValid) - fix(fs / 10)):numel(idxValid);
        idxValid(idxPre) = 0;
        idxValid(idxPost) = 0;
        idxValid = logical(idxValid);

        % filter original state
        currSBeta1 = abs(hilbert(filtfilt(bpFilt_Beta1, currS')))';
        currSBeta2 = abs(hilbert(filtfilt(bpFilt_Beta2, currS')))';
        currSBeta = (currSBeta1 + currSBeta2) / 2;
        
        % filter next state
        currS_stimBeta1 = abs(hilbert(filtfilt(bpFilt_Beta1, currS_stim')))';
        currS_stimBeta2 = abs(hilbert(filtfilt(bpFilt_Beta2, currS_stim')))';
        currS_stimBeta = (currS_stimBeta1 + currS_stimBeta2) / 2;
    
        % reward calculation
        % if based on power only
        if contains(strRewardMethod, 'power')
            currReward = sum(currSBeta(idxValid) - currS_stimBeta(idxValid));
        
        % if based on beta burst
        elseif contains(strRewardMethod, 'burst')
            burstDurationS = sum(currSBeta(idxValid) > currBurstThresh);
            burstDurationS_stim = sum(currS_stimBeta(idxValid) > currBurstThresh);
            
            currReward = burstDurationS - burstDurationS_stim;
        end

        % append to outer list
        reward_cond{i} = currReward;
    end

    % convert back to array
    reward_cond = cat(2, reward_cond{:});

    % if wish to normalize
    if contains(strRewardMethod, 'Norm')
        reward_cond = reward_cond / std(reward_cond);
    end

    % append to outer structure
    C{condsel} = reward_cond;

end

% save reward as output
if boolSave
    if init
        % initial save
        outputPath = fullfile(R.DRL.outPath, 'initBuffer');
        if ~exist(outputPath, 'dir')
            mkdir(outputPath);
        end

        % save output structure
        strOutputName = sprintf('reward_%s', strRewardMethod);
        save(fullfile(outputPath, strOutputName), 'C', '-v7.3');
    else
        % part of training
        outputPath = fullfile(R.DRL.outPath, strSubFolder');

        % save output structure
        strOutputName = sprintf('reward_%s', strRewardMethod);
        save(fullfile(outputPath, strOutputName), 'C', '-v7.3');
    end
end

end
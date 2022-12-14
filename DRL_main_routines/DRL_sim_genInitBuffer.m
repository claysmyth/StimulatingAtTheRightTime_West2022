function [X, X_stim, A, R] = DRL_sim_genInitBuffer(R, varargin)
%% Input parsing
% Handle the optional inputs
p = inputParser;
p.KeepUnmatched = true;

addParameter(p, 'SScomb', 1, ...                                            % stimulation type
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'stimAmp', 1e3, ...                                         % stimulation amplitude
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

parse(p,varargin{:});

% Handles incorrect inputs
UnmatchedParam = fieldnames(p.Unmatched);
if ~isempty(UnmatchedParam)
    error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
end

% unpacking variable
SScomb = p.Results.SScomb;
stimAmp = p.Results.stimAmp;

%% Model-based rollout
% Now load model and run simulation without external stim
% generate the state variable
% tic;
fprintf('\n>>Simulating the original states\n')
[R, m, p] = loadModelForDRL(R);                                             % load saved model fitted on rat data

u = simNoiseIntr(R, m);                                                     % obtain the initial intrinsic noise
[~, ~, ~, xsims, xsims_gl, ~, R] = computeSimData120319(R, m, u, p, 0);
% timeSim = toc;

%% Perform epoching

% update parameter set R
% tic;
R = updateREpochForDRL(R);

% now perform epoching in state and random noise variable
fprintf('\n>>Epoching data\n')
X = DRL_epochData(xsims, xsims_gl, u, R);
% timeEpoch = toc;

%% next test with simulation

% run simulation with fixed seed to figure out eps
% tic;
fprintf('\n>>Estimating threshold for sigmoid\n')
[eps, ~] = estimateEps(R);
R.IntP.phaseStim.eps = eps;

% now run simulation with actual stim
% probably linspace between 1e3 and 2e3
fprintf('\n>>Simulating random action\n')
[X_stim, A, R] = DRL_resim_bufferStim(X, R, m, p, ...
    'SScomb', SScomb, 'stimAmp', stimAmp);

% also downsample all data
[X, X_stim, A] = DRL_downSampleData(X, ...
    X_stim, A, R);
% timeResim = toc;

end
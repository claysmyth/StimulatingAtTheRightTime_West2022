function [X, R] = DRL_sim_genInitEvalBuffer(R)
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

%% Downsample data (redacted version)

X = DRL_downSampleData_Redact(X, R);

end
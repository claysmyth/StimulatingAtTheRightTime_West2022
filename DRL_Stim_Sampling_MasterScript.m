%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

fprintf('\n>>Initializing models\n')
R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

%% Model-based rollout
% Now load model and run simulation without external stim
% generate the state variable
tic;
fprintf('\n>>Simulating the original states\n')
[R, m, p] = loadModelForDRL(R);                                             % load saved model fitted on rat data
u = simNoiseIntr(R, m);                                                     % obtain the initial intrinsic noise
[~, ~, ~, xsims, xsims_gl, wflag, R] = ...                                     % generate the state variable s
    computeSimData120319(R,m,u,p,0);
timeSim = toc;

%% Perform epoching

% update parameter set R
tic;
R = updateREpochForDRL(R);

% now perform epoching in state and random noise variable
fprintf('\n>>Epoching data\n')
X = DRL_epochData(xsims, xsims_gl, u, R);
timeEpoch = toc;

%% next test with simulation

% run simulation with fixed seed to figure out eps
fprintf('\n>>Estimating threshold for sigmoid\n')
[eps, ~] = estimateEps(R);
R.IntP.phaseStim.eps = eps;

% now run simulation with actual stim
% probably linspace between 1e3 and 2e3
fprintf('\n>>Simulating random action\n')
[X_stim, A, R] = DRL_sim_bufferStim(X, R, m, p, 'SScomb', 2);

[~, ~, ~] = DRL_sim_bufferStim(X, R, m, p, 'SScomb', 3, 'exStimStruct', A);
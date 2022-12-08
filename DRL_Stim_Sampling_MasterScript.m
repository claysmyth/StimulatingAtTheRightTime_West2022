%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

%% Model-based rollout
% Now load model and run simulation without external stim
% generate the state variable
tic;
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
X = DRL_epochData(xsims, xsims_gl, u, R);
timeEpoch = toc;

%% next test with simulation




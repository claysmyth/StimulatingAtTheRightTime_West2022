%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

%% Now load model and run simulation without external stim

[R, m, p] = loadModelForDRL(R);                                             % load saved model fitted on rat data
u = simNoiseIntr(R, m);
[~, pnew, feat_sim, xsims, xsims_gl, wflag] = ...
    computeSimData120319(R,m,u,p,0);

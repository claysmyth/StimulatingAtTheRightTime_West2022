% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%      Master Script for Simulations/Figures in West et al. (2022)      %%
%       "Stimulating at the Right Time to Recover Network States          %
%          in a model of the Cortico-Basal Ganglia-Thalamic Circuit"      %
%                                                                         %
% This is the main script that reproduces the figures in the              %
% manuscript. This analysis takes forward a model of the cortico-basal    %
% ganglia circuit implemented in Dynamic Causal Modelling and previously  %
% described in van Wijk et al. (2018). This model was fit to data from    %
% West al. (2018) using an optimization based upon Approximate Bayesian   %
% Computation and described in West et al. (2021). Many of these analyses %
% are based upon scripts included in publicly available toolboxes and     %
% generously provided by their respective authors. Please see             %
% 'ABC_dependencies' for their respective licenses.                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Timothy West; MRC Brain Networks Dynamics Unit,Nuffield Department of   %
% Clinical Neurosciences, University of Oxford;                           %
% Wellcome Centre for Human Neuroscience, University College London.      %
% 2018-2022                                                               %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% clean up environment
restoredefaultpath
clear; close all
% set plot defaults
set(0,'defaultAxesFontSize',16)

%% Initialise Paths and Configuration Structure (R)
R = project_AddPaths();             % Organise file paths and add dependencies to path
R = setupBasalGangliaModel(R);      % Add configurations/settings to R

%% Figure (1): Model fit and comparison with data
BAA_plotModelFit(R);                % This plots the model

%% Figure (2): Sweep over connections and plot spectra
BAA_sim_ConnectionSweep(R)       % Perform simulations with connectivity Sweep
BAA_plotSweepSpectraBasic(R);       % Plots the spectra from the sweep

%% Figure (3): Define and Analyse Bursts in Spontaneous Data
BAA_computeBurstWrapper(R);      % Computes bursts and saves spontaneous network states
BAA_computeBurstWrapper_emp(R); % Computes bursts  and saves spontaneous network states for animal data

%% Figure (4): Model of Stimulation
BAA_sim_phaseLockedStim(R);         % Do the simulations       
BAA_computeStimAnalysis_sweep(R,1); % Sweep across conditions and compute features
BAA_plotStimTracesBase(R);          % Plot stim outcomes in burst model

%% Figure(5): Compare the states
BAA_fingerprintCompare_Sweep(R);    % Fingerprint comparison between stim and spontaneous
BAA_fingerprintCompare_emp(R);      % Fingerprint comparison between stim and spontaneous

%% Supplementary Information
% Investigate the role of Obs Noise on Phase Estimation
BAA_sim_phaseLockedStim_SICompObsNoise(R); % this has done SScomb =1

% Plot different stim types
BAA_plotStimTypeOutcomes(R); % SI- Plot stim types

% Threshold and Sim
BAA_sim_phaseLockedStim_SICompEps(R); % Effect of threshold on stim effects
BAA_plotStimTypeOutcomes_Eps(R)

% Validation of phase estimation technique - zero crossing
peValidation
peValidation_sweepObservationNoise % this one looks at observation noise

%% These functions need the raw experimental data
% % plotDataComparison % plots the time series in figure 1A and C
% % % Bimodal data inspection
% % prepareRatData_3Gauss_Group_NPD(R)
% % prepareRatData_3Gauss_Subject_NPD(R)

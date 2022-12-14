function [eps, R] = estimateEps(Rorg)
% estimate the eeps (threshold) value from loaded model

%% Load in the model
load(fullfile(Rorg.rootn, 'data', 'modelfit', 'SimModelData_M10.mat'), ...
    'R','m','p');
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
porg = p;

% fake stim parameters
% Stimulating  STN - STN cDBS- DC modulation
% R.IntP.phaseStim.stimFx = @dcStim_v1;
senssite = 4; % M2
stimsite = 1; % STN
R = typeIstimPars_v3(R);
% R.IntP.phaseStim.stimAmp = 1e2; % 1e4;%
% R.IntP.phaseStim.upperiod = 15;%
% R.IntP.phaseStim.stimlength = 15;%

% Setup stim parameters
R.IntP.phaseStim.sensStm = [senssite stimsite];

%% settting model params

% set the simulation hyperparameters
R.frqz = 6:.2:148;
% Simulation Conditions
R = setSimTime(R,32); %128);
    
% Trans Options
R.obs.trans.norm = 0;
R.obs.gainmeth = {};

% Give all timeseries the same input - makes comparable
rng(6543)
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
XBase = porg;
    
% integration function to be used
R.IntP.intFx = @spm_fx_compile_120319_stim;


%% Setup Base Model
Pbase = XBase;

% Simulate Base Model
uc_ip{1} = uc;
R.IntP.phaseStim.switch = 0 ;
R.IntP.phaseStim.phaseshift = 0;
R.IntP.compFx = @nullComp;
[~,~,~,xsim_gl,~] = ...
    computeSimData120319(R,m,uc_ip{1},Pbase,0, 0, 0);

% Work out the threshold
R.IntP.phaseStim.eps = 0;
[~,R] = zeroCrossingPhaseStim_v3([],R,0,xsim_gl{1},R.IntP.dt);
eps = R.IntP.phaseStim.eps;

end
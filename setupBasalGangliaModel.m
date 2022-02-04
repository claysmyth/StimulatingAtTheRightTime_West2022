function R = setupBasalGangliaModel(R)
% Set plotting defaults
set(0,'defaultAxesFontName','Arial')

% File managment
R.projectn = 'Rat_NPD';
R.out.tag = 'InDrt_ModCompRev2';
R.filepathn = [R.rootn 'data\storage'];

% 
%% DATA SPECIFICATION
R.data.datatype = 'NPD'; %%'NPD'
R.frqz = [6:.2:68];
R.frqzfull = [1:.2:200]; % used for filters
R.chloc_name = {'MMC','STR','GPE','STN'};
R.chsim_name = {'MMC','STR','GPE','STN','GPI','THAL'};
R.condnames = {'OFF'};
% Spectral characteristics
R.obs.csd.df = 0.5;
R.obs.csd.reps = 32;

%% INTEGRATION
% Main dynamics function
R.IntP.intFx = @spm_fx_compile_120319;
R.IntP.compFx= @compareData_100717;

R.IntP.dt = .0005;
R.IntP.Utype = 'white_covar';
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays

N = R.obs.csd.reps; % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df))); % order of NPD for simulated data
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;
R.IntP.nt = R.IntP.tend/R.IntP.dt;
R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);

dfact = fsamp/(2*2^(R.obs.SimOrd));
disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));

%% OBSERVATION 
% observation function
R.obs.obsFx = @observe_data;
R.obs.gainmeth = {'unitvar'};
R.obs.brn =3; % % burn in time
LF = [1 1 1 1]*10; % Leadfields
R.obs.LF = LF;

% Data Features
% fx to construct data features
R.obs.transFx = @constructNPDMat_190618;
% These are options for transformation (NPD)
R.obs.glist =0;
R.obs.trans.logdetrend =0;
R.obs.trans.norm = 0;
R.obs.trans.gauss = 0;
R.obs.logscale = 0;


%% OBJECTIVE FUNCTION
R.objfx.feattype = 'ForRev'; % Asymetric
R.objfx.specspec = 'cross';  % include cross terms

%% PLOTTING
R.plot.outFeatFx = @npdplotter_110717; %%@;csdplotter_220517
R.plot.save = 'False';


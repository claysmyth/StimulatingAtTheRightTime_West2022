function R = setupBasalGangliaModel(R)
% Set plotting defaults
set(0,'defaultAxesFontName','Arial')

% File managment
R.projectn = 'DRL_Stim';                                                    % set project name
R.out.tag = 'StimSampleBuffer_v1';                                          % set output buffer
R.filepathn = fullfile(R.rootn, 'data', 'storage');                         % set output path

% 
%% DATA SPECIFICATION
R.data.datatype = 'NPD'; %%'NPD'
R.frqz = [6:.2:68];                                                         % frequency range of interest
R.frqzfull = [1:.2:200]; % used for filters                                 % full frequency range of the data
R.chloc_name = {'MMC','STR','GPE','STN'};
R.chsim_name = {'MMC','STR','GPE','STN','GPI','THAL'};                      % name of simulated channels
R.condnames = {'OFF'};                                                      % name of conditions

% Spectral characteristics
R.obs.csd.df = 0.5;                                                         % degree of freedom for CSD, Dont Care
R.obs.csd.reps = 32;                                                        % number of epochs for freq resolution, Dont Care

%% INTEGRATION
% Main dynamics function
R.IntP.intFx = @spm_fx_compile_120319;                                      % integration function for simulation
%R.IntP.intFx = @spm_fx_compile_DRL;  
R.IntP.compFx= @compareData_100717;                                         

R.IntP.dt = .0005;                                                           % current sampling rate at 1kHz
R.IntP.Utype = 'white_covar';                                               % sets the background noise generation method
R.IntP.buffer = ceil(0.050*(1/R.IntP.dt)); % buffer for delays              % length of buffer of length 50ms for delays

N = R.obs.csd.reps;                                                         % Number of epochs of desired frequency res
fsamp = 1/R.IntP.dt;                                                        % compute the sampling rate
R.obs.SimOrd = floor(log2(fsamp/(2*R.obs.csd.df)));                         % order of NPD for simulated data, Dont Care
R.IntP.tend = (N*(2^(R.obs.SimOrd)))/fsamp;                                 % length of recording in seconds             
R.IntP.nt = R.IntP.tend/R.IntP.dt;                                          % number of time samples in the recording          
R.IntP.tvec = 0:R.IntP.dt:R.IntP.tend;                            % time stamp vector

dfact = fsamp/(2*2^(R.obs.SimOrd));                                         % actual stimulation df (freq res)
disp(sprintf('The target simulation df is %.2f Hz',R.obs.csd.df));
disp(sprintf('The actual simulation df is %.2f Hz',dfact));

%% OBSERVATION 
% observation function
R.obs.obsFx = @observe_data;                                                % function handle for getting LFP from NMM model
R.obs.gainmeth = {'unitvar'};                                               % method for controlling channel gain
R.obs.brn =3; %                                                             % burn in time in seconds
LF = [1 1 1 1]*10;                                                          % Leadfields
R.obs.LF = LF;

% Data Features
% fx to construct data features
R.obs.transFx = @constructNPDMat_190618;                                    % transformation function for NPD, Dont Care
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

%% DRL STIM PARAMS

R.DRL.epLen = 1;                                                            % length of each sample/epoch in seconds
R.DRL.nEpLen = R.DRL.epLen / R.IntP.dt;                                     % number of samples in each sample/epoch
R.DRL.dsFs = 200;                                                           % desired sampling rate of downsampled signal in Hz

% set output path
U = char(java.lang.System.getProperty('user.name'));
switch U
    case 'jyao'
        R.DRL.outPath = fullfile('/home', 'jyao', 'local', 'data', 'starrlab', 'DRL');
        if ~exist(R.DRL.outPath, 'dir')
            mkdir(R.DRL.outPath);
        end

    otherwise
        error("need ot specify output path for saving data")

end

% sanity check
assert(isfield(R.DRL, 'outPath'));

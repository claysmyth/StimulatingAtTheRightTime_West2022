function [R] = DRL_sim_bufferStim(Rorg, varargin)
close all
%% Load in the model
load(fullfile(Rorg.rootn, 'data', 'modelfit', 'SimModelData_M10.mat'), ...
    'R','m','p');
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
porg = p;

%stimlength = 15 % TODO: Determine if this shoud be an input parameter?

for SScomb = 1 % These are the different stim types
    %% Define stimulation conditions
    if SScomb == 1
        % Stimulating  STN - STN cDBS- DC modulation
        senssite = 1; % M2
        stimsite = 4; % STN
        stimFx = @dcStim_v1;
        % phflag = 0;
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = 1e2; % 1e4;%
        R.IntP.phaseStim.upperiod = 15;%
        %R.IntP.phaseStim.stimlength = stimlength;%
    elseif SScomb == 2
        % Stimulating  STN - STN DBS with random DC modulation
        senssite = 1; % M2
        stimsite = 4; % STN
        stimFx = @getRandomStimulationAction;
        % phflag = 0;
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = 1e2; % 1e4;%
        R.IntP.phaseStim.upperiod = 15;%
        %R.IntP.phaseStim.stimlength = stimlength;%
    elseif SScomb == 3
        % Stimulating STN - STN DBS with externally (e.g. DRL output) determined stim vector
        senssite = 1; % M2
        stimsite = 4; % STN
        stimFx = @passExternalStimVector % External stim vector
        R.IntP.phaseStim.stimAmp = 1e2; % 1e4;%
        R.IntP.phaseStim.upperiod = 15;%
        %R.IntP.phaseStim.stimlength = stimlength;% 
    end
end

R.IntP.phaseStim.sensStm = [senssite stimsite];

%% Connection Sets
rootan = fullfile(Rorg.rootn, 'data', 'ConnectionSweep');
load(fullfile(rootan, ['BB_' Rorg.out.tag '_ConnectionSweep_CON_1_ck_1.mat']), ...
    'ck_1');                                                            % load connection bands (CON 1 and 2 have same data)

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

% Phase To be Tested
R.IntP.intFx = @spm_fx_compile_120319_stim;

% Running only DC cDBS paradigms. No phase information.
phaseShift = 0;

%% Loop through Connections
for CON =1:2
    feat_sim_save = {};
    xsim_ip = {};
    for state = 1:size(ck_1,2)
        %% Setup Base Model
        Pbase = XBase;
        if CON == 1 % Hyperdirect
            Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(CON,state)); %
        elseif CON == 2 % Pallidal-subthalamo
            Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_1(CON,state)); %
        end
        
        % Simulate Base Model
        uc_ip{1} = uc;
        R.IntP.phaseStim.switch = 0 ;
        R.IntP.phaseStim.phaseshift = 0;
        R.IntP.compFx = @nullComp;
        [~,~,feat_sim_base{1},xsim_gl,xsim_ip_base{1}] = computeSimData120319(R,m,uc_ip{1},Pbase,0);
        
        % Work out the threshold
        R.IntP.phaseStim.eps = 0;
        [~,R] = zeroCrossingPhaseStim_v3([],R,0,xsim_gl{1},R.IntP.dt);
        eps(state) =  R.IntP.phaseStim.eps;
        
        %% Do Stimulation
        R.IntP.phaseStim.switch = 1;
        xsim_ip_stim = cell(1,12); feat_sim_stim = cell(1,12); pU = cell(1,12);
        
        Rpar = R;

        % Modulate the phase (setting to 0 in this case)
        Rpar.IntP.phaseStim.phaseshift = phaseShift;
        
        Rpar.IntP.phaseStim.stimFx = stimFx;
        % Simulate with Stimulation
        [~,~,feat_sim_stim{p},~,xsim_ip_stim{p},~]  = computeSimData120319(Rpar,m,uc_ip{1},Pbase,0);

        uexs = load([Rpar.rootn 'data/phaseStimSave/stim_tmp_' sprintf('%3.f',1000*Rpar.IntP.phaseStim.phaseshift)],'uexs');
        pU{p} = uexs.uexs(Rpar.IntP.phaseStim.sensStm(2),round(Rpar.obs.brn*(1/R.IntP.dt))+1:end);
        disp([CON state p])
        
        rmdir([R.rootn 'data/phaseStimSave/'],'s')
        % Create save version
        feat_sim_save{1,state} = feat_sim_base; feat_sim_save{2,state} = feat_sim_stim;
        xsim_ip{1,state} = xsim_ip_base; xsim_ip{2,state} = xsim_ip_stim;
        pU_save{state} = pU;
        
    end
    %         Rorg.out.tag = 'tester';

    % TODO: DECIDE HOW TO SAVE NEW STIM+SIM DATA
    rootan = [Rorg.rootn 'data/<FileNameHere>'];
    mkdir(rootan)
    
    save([rootan '/DRL_' Rorg.out.tag '_rerun_sim_with_stim_' num2str(CON) '_feat' num2str(SScomb) '.mat'],'feat_sim_save')
    save([rootan '/DRL_' Rorg.out.tag '_rerun_sim_with_stim_' num2str(CON) '_xsim' num2str(SScomb) '.mat'],'xsim_ip','-v7.3')
    save([rootan '/DRL_' Rorg.out.tag '_rerun_sim_with_stim_' num2str(CON) '_Rout' num2str(SScomb) '.mat'],'R')
    save([rootan '/DRL_' Rorg.out.tag '_rerun_sim_with_stim_' num2str(CON) '_pU_save' num2str(SScomb) '.mat'],'pU_save')
end
disp(SScomb)
function [R] = DRL_sim_bufferStim(X, Rorg, morg, porg, varargin)
% In current iteration, varargin is an optional parameter for passing in 
% an external stimulation vector. Size(varargin) = 1.

%% Input parsing

% Handle the optional inputs
p = inputParser;
p.KeepUnmatched = true;

addParameter(p, 'external_stim_vector', {});                                % optional handle for passing in external stim

addParameter(p, 'stimlength', 15, ...                                       % length of stimulation
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'upperiod', 15, ...                                       % length of stimulation
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'SScomb', 1, ...                                            % stimulation type
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'stimAmp', 1e2, ...                                         % stimulation amplitude
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

parse(p,varargin{:});

% Handles incorrect inputs
UnmatchedParam = fieldnames(p.Unmatched);
if ~isempty(UnmatchedParam)
    error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
end

% unpacking variable
external_stim_vector = p.Results.external_stim_vector;
stimlength = p.Results.stimlength;
upperiod = p.Results.upperiod;
SScomb = p.Results.SScomb;
stimAmp = p.Results.stimAmp;

%% Load in the model

R = deepCopyStruct(Rorg, 1);

for i = 1 % These are the different stim types
    %% Define stimulation conditions
    if SScomb == 1
        % Stimulating  STN - STN cDBS- DC modulation
        senssite = 1; % M2
        stimsite = 4; % STN
        R.IntP.phaseStim.stimFx = @dcStim_v1;
        % phflag = 0;
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = stimAmp; % 1e4;%
        R.IntP.phaseStim.upperiod = upperiod;%
        R.IntP.phaseStim.stimlength = stimlength;%
    
    elseif SScomb == 2
        % Stimulating  STN - STN DBS with random DC modulation
        senssite = 1; % M2
        stimsite = 4; % STN
        R.IntP.phaseStim.stimFx = @getRandomStimulationAction;
        % phflag = 0;
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = stimAmp; % 1e4;%
        R.IntP.phaseStim.upperiod = upperiod;%
        R.IntP.phaseStim.stimlength = stimlength;%
    
    elseif SScomb == 3
        % Stimulating STN - STN DBS with externally (e.g. DRL output) 
        % determined stim vector
        senssite = 1; % M2
        stimsite = 4; % STN
        R.IntP.phaseStim.stimFx = @passExternalStimVector;               % External stim vector
        
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = stimAmp; % 1e4;%
        R.IntP.phaseStim.upperiod = upperiod;%
        R.IntP.stimtype = 'external';
        R.IntP.phaseStim.stimlength = stimlength;% 
        R.IntP.phaseStim.externalStim = external_stim_vector;
    end
end

R.IntP.phaseStim.sensStm = [senssite stimsite];
Rorg = deepCopyStruct(R, 1);

%% Reset some flags for simulation

% reshuffle random seed
rng shuffle;

% set integration function to simulation with external stim
R.IntP.intFx = @spm_fx_compile_DRL_stim;

% set some spectra flags, Don't care
R.obs.csd.df = 0.5;
R.frqz = 6:.2:148;

% switch on stimulation
R.IntP.phaseStim.switch = 1;
R.IntP.phaseStim.phaseshift = 0;                                            % Running only DC cDBS paradigms. No phase information.

%% Do stimulation

% loop through the different conditions interested in
for condsel = 1:numel(R.condnames)
    % unpack X variable
    X_cond = X{condsel};
    numEpoch = X_cond.metadata.numEpoch;

    % next generate numEpoch copies of R, m and p
    vecR = deepCopyStruct(R, numEpoch);
    vec_m = deepCopyStruct(m, numEpoch);
    vec_p = deepCopyStruct(p, numEpoch);

    % further unpack for parfor
    for i = 1:numEpoch
        vecUFull{i} = {[X_cond.u.Buffer{i}; X_cond.u.SPrev{i}; ...
            X_cond.u.S{i}]};
        tvecRel_temp = [X_cond.xsims_gl.tvec_Buffer{i}, ...
            X_cond.xsims_gl.tvec_SPrev{i}, X_cond.xsims_gl.tvec_S{i}];
        vec_tvecRel{i} = tvecRel_temp - X_cond.xsims_gl.tvec_S{i}(1);
    end
    xsimsBuff = X_cond.xsims.Buffer;
    
    % now initialize the parallelized for loop
    for i = 1:numEpoch
        % unpack meta structure
        R_i = vecR{i};
        m_i = vec_m{i};
        p_i = vec_p{i};

        % form arrays necessary for resimulation
        uFull_i = vecUFull{i};                                              % full intrinsic noise
        xsimsBuff_i = xsimsBuff{i};                                         % initial states for simulation
        tStepEnd_i = size(uFull_i{1}, 1);                                   % end time step for simulation
        tvecRel_i = vec_tvecRel{i};
        
        % re-run simulation
        [xsims_stim, xsim_gl_stim] = ...
            computeSimDataDRL_wUOnly(R_i, xsimsBuff_i, ...
            m_i, uFull_i, p_i, tStepEnd_i, tvecRel_i);

        % append to outer array
        t1 = 1;


    end

end

% obtain number of epochs and make copy of hyper structure


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
        [~,~,feat_sim_base{1},xsim_gl,xsim_ip_base{1}] = computeSimDataDRL(R,m,uc_ip{1},Pbase,0, external_stim_vector);
        
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
        [~,~,feat_sim_stim{p},~,xsim_ip_stim{p},~]  = computeSimDataDRL(Rpar,m,uc_ip{1},Pbase,0, external_stim_vector);

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
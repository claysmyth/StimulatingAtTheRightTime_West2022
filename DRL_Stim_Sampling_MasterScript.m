%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

%% Model-based rollout
% Now load model and run simulation without external stim
% generate the state variable
[R, m, p] = loadModelForDRL(R);                                             % load saved model fitted on rat data
u = simNoiseIntr(R, m);                                                     % obtain the initial intrinsic noise
[~, ~, ~, xsims, xsims_gl, wflag, R] = ...                                     % generate the state variable s
    computeSimData120319(R,m,u,p,0);
R.obs.obsFx = @DRL_observe_data;

% now perform epoching in state and random noise variable
% first determine the number of epochs based on brn time and
% reset sim time
R.DRL.brnNEpoch = ceil(R.obs.brn / R.DRL.epLen);                            % number of epochs per burn in time
R.DRL.brnNSample = R.obs.brn / R.IntP.dt;                                   % number of samples in epochs to discard
R.DRL.tendReCalc = floor(R.IntP.tend / R.DRL.epLen) * R.DRL.epLen; 
if R.DRL.tendReCalc ~= R.IntP.tend                         
    R.IntP.tend = R.DRL.tendReCalc;                                             
    R.IntP.nt = R.IntP.tend/R.IntP.dt;      
    R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);
    R.IntP.tvec_obs(R.IntP.tvec_obs > R.IntP.tend) = [];
end

% next obtain the time vectors from the simulated R vector
tvec = R.IntP.tvec;                                                         % time for xsim, all state vectors [0:dt:tend]
tvecObs = R.IntP.tvec_obs;                                                  % time for xsim_gl, observed LFP [0, brn:dt:tend]
tvecU = tvec(1:(size(tvec, 2) - 1));

% loop through all condition names
for condsel = 1:numel(R.condnames)
    assert(size(xsims{condsel}, 2) == size(tvec, 2))
    assert(size(xsims_gl{condsel}, 2) == size(tvecObs, 2))
    assert(size(u{condsel}, 1) == size(tvecU, 2))

    % perform epoching in xsims, xsims_gl and u vectors
    % start with xsims, note here add in one extra second for
    % loading in right buffer for initializing all transitions
    tStart_xsims_Full = (R.obs.brn + 1):R.DRL.epLen:R.IntP.tend;
    idxEpochCount = 1;
    for i = (R.DRL.brnNEpoch + 1):(numel(tStart_xsims_Full) - 1)
        % obtain the start and end time for now
        tStartCurr = tStart_xsims_Full(i);
        tEndCurr = tStart_xsims_Full(i + 1);
        tStartPrevCurr = tStart_xsims_Full(i - R.DRL.brnNEpoch);

        % perform epoching in state
        idxSCurr = tvec >= tStartCurr & tvec < tEndCurr;
        xsims_Cond_S{idxEpochCount} = xsims{condsel}(:, idxSCurr);
        tvec_Cond_S{idxEpochCount} = tvec(idxSCurr);

        % obtain all previous states (for settling behavior)
        idxSPrevCurr = tvec >= tStartPrevCurr & tvec < tStartCurr;
        xsims_Cond_SPrev{idxEpochCount} = xsims{condsel}(:, idxSPrevCurr);
        tvec_Cond_SPrev{idxEpochCount} = tvec(idxSPrevCurr);

        % obtain buffer of data for starting simulation
        idxStartPrevCurr = find(idxSPrevCurr, 1, 'first');
        idxStartBufferCurr = idxStartPrevCurr - R.IntP.buffer;
        xsims_Cond_Buffer{idxEpochCount} = ...
            xsims{condsel}(:, idxStartBufferCurr:(idxStartBufferCurr+R.IntP.buffer-1));
        tvec_Cond_Buffer{idxEpochCount} = ...
            tvec(idxStartBufferCurr:(idxStartBufferCurr+R.IntP.buffer-1));

        % advance index of epoch
        idxEpochCount = idxEpochCount + 1;
    end

    % next perform epoching in xsim_gl
    tStart_xsims_gl_Full = (R.obs.brn + 1):R.DRL.epLen:R.IntP.tend;
    idxEpochCount = 1;
    for i = (R.DRL.brnNEpoch + 1):(numel(tStart_xsims_gl_Full) - 1)
        % obtain the start and end time for now
        tStartCurr = tStart_xsims_gl_Full(i);
        tEndCurr = tStart_xsims_gl_Full(i + 1);
        tStartPrevCurr = tStart_xsims_gl_Full(i - R.DRL.brnNEpoch);
    
        % perform epoching in state
        idxSObsCurr = tvecObs >= tStartCurr & tvecObs < tEndCurr;
        xsims_gl_Cond_S{idxEpochCount} = xsims_gl{condsel}(:, idxSObsCurr);
        tvecObs_Cond_S{idxEpochCount} = tvecObs(idxSObsCurr);

        % obtain all previous states (for settling behavior)
        idxSObsPrevCurr = tvecObs >= tStartPrevCurr & tvecObs < tStartCurr;
        xsims_gl_Cond_SPrev{idxEpochCount} = xsims_gl{condsel}(:, idxSObsPrevCurr);
        tvecObs_Cond_SPrev{idxEpochCount} = tvecObs(idxSObsPrevCurr);

        % obtain buffer of data for comparison purposes
        idxStartPrevCurr = find(idxSObsPrevCurr, 1, 'first');
        idxStartBufferCurr = idxStartPrevCurr - R.IntP.buffer;
        xsims_gl_Cond_Buffer{idxEpochCount} = ...
            xsims_gl{condsel}(:, idxStartBufferCurr:(idxStartBufferCurr+R.IntP.buffer-1));
        tvecObs_Cond_Buffer{idxEpochCount} = ...
            tvecObs(idxStartBufferCurr:(idxStartBufferCurr+R.IntP.buffer-1));


        % advance index of epoch
        idxEpochCount = idxEpochCount + 1;
    end

    % finally perform epoching in u - intrinsic noise
    tStart_u_Full = (R.obs.brn + 1):R.DRL.epLen:R.IntP.tend;
    idxEpochCount = 1;
    for i = (R.DRL.brnNEpoch + 1):(numel(tStart_u_Full) - 1)
        % obtain the start and end time for now
        tStartCurr = tStart_u_Full(i);
        tEndCurr = tStart_u_Full(i + 1);
        tStartPrevCurr = tStart_u_Full(i - R.DRL.brnNEpoch);

        % perform epoching in state
        idxUCurr = tvecU >= tStartCurr & tvecU < tEndCurr;
        u_Cond_S{idxEpochCount} = u{condsel}(idxUCurr, :);
        tvecU_Cond_S{idxEpochCount} = tvecU(idxUCurr);

        % obtain all previous states (for settling behavior)
        idxUPrevCurr = tvecU >= tStartPrevCurr & tvecU < tStartCurr;
        u_Cond_SPrev{idxEpochCount} = u{condsel}(idxUPrevCurr, :);
        tvecU_Cond_SPrev{idxEpochCount} = tvecU(idxUPrevCurr);

        % obtain buffer of data for starting simulation
        idxStartPrevCurr = find(idxUPrevCurr, 1, 'first');
        idxStartBufferCurr = idxStartPrevCurr - R.IntP.buffer;
        u_Cond_Buffer{idxEpochCount} = ...
            u{condsel}(idxStartBufferCurr:(idxStartBufferCurr+R.IntP.buffer-1), :);
        tvecU_Cond_Buffer{idxEpochCount} = ...
            tvecU(idxStartBufferCurr:(idxStartBufferCurr+R.IntP.buffer-1));

        % advance index of epoch
        idxEpochCount = idxEpochCount + 1;
    end
end

% next test with simulation
% make deep copy of R struct
fnR = fieldnames(R);
for j = 1:numel(fnR)
    R_i.(fnR{j}) = R.(fnR{j});
end

% make deep copy of m struct
fnM = fieldnames(m);
for j = 1:numel(fnM)
    m_i.(fnM{j}) = m.(fnM{j});
end

fnP = fieldnames(p);
for j = 1:numel(fnP)
    p_i.(fnP{j}) = p.(fnP{j});
end

% piece together the u variable
i = 1;
uFull = {[u_Cond_Buffer{i}; u_Cond_SPrev{i}; u_Cond_S{i}]};
xBuff = xsims_Cond_Buffer{i};
xComp = [xsims_Cond_SPrev{i}, xsims_Cond_S{i}];
xCompFull = [xBuff, xComp];

xGLCompFull = [xsims_gl_Cond_Buffer{i}, xsims_gl_Cond_SPrev{i}, xsims_gl_Cond_S{i}];
tStepEnd = size(uFull{1}, 1);
tvecRel_i = [tvecObs_Cond_Buffer{i}, tvecObs_Cond_SPrev{i}, tvecObs_Cond_S{i}];
tvecRel_i = tvecRel_i - tvecObs_Cond_S{i}(1);

% now engineer the R struct
[xsimsAlt dum wflag] = spm_fx_compile_120319_debug(R_i, xBuff, uFull, p_i, m_i, xComp, tStepEnd);

% Run Observer function
% Subloop is local optimization of the observer gain
glorg = p_i.obs.LF;
gainlist = R_i.obs.glist;
for gl = 1:length(gainlist)
    % run observation function
    p_i.obs.LF = glorg+gainlist(gl);
    if isfield(R_i.obs,'obsFx')
        [xsims_glAlt{gl}, ~, wflag(1)] = R_i.obs.obsFx(xsimsAlt,m_i,p_i,R_i, tvecRel_i);
    else
        xsims_glAlt{gl} =xsimsAlt;
    end
end
xsims_glAlt = xsims_glAlt{1};

t1 = 1;

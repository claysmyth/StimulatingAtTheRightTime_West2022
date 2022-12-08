function X = DRL_epochData(xsims, xsims_gl, u, R)

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
    %% epoch xsims

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
    
    %% epoch xsims_gl

    % next perform epoching in xsims_gl
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
    
    %% epoch U

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

    %%  append to outer structure;

    % create cond-specific struct
    % first append all state-related variables
    X_Cond.xsims.Buffer = xsims_Cond_Buffer;
    X_Cond.xsims.SPrev = xsims_Cond_SPrev;
    X_Cond.xsims.S = xsims_Cond_S;

    X_Cond.xsims.tvec_Buffer = tvec_Cond_Buffer;
    X_Cond.xsims.tvec_SPrev = tvec_Cond_SPrev;
    X_Cond.xsims.tvec_S = tvec_Cond_S;
    X_Cond.xsims.tStart = tStart_xsims_Full;

    % next append all LFP-related variables
    X_Cond.xsims_gl.Buffer = xsims_gl_Cond_Buffer;
    X_Cond.xsims_gl.SPrev = xsims_gl_Cond_SPrev;
    X_Cond.xsims_gl.S = xsims_gl_Cond_S;

    X_Cond.xsims_gl.tvec_Buffer = tvecObs_Cond_Buffer;
    X_Cond.xsims_gl.tvec_SPrev = tvecObs_Cond_SPrev;
    X_Cond.xsims_gl.tvec_S = tvecObs_Cond_S;
    X_Cond.xsims_gl.tStart = tStart_xsims_gl_Full;

    % finally append all u-related variables
    X_Cond.u.Buffer = u_Cond_Buffer;
    X_Cond.u.SPrev = u_Cond_SPrev;
    X_Cond.u.S = u_Cond_S;

    X_Cond.u.tvec_Buffer = tvecU_Cond_Buffer;
    X_Cond.u.tvec_SPrev = tvecU_Cond_SPrev;
    X_Cond.u.tvec_S = tvecU_Cond_S;
    X_Cond.u.tStart = tStart_u_Full;

    % also append metadata (for debugging)
    X_Cond.metadata.tvec = tvec;
    X_Cond.metadata.tvecObs = tvecObs;
    X_Cond.metadata.tvecU = tvecU;
    X_Cond.metadata.numEpoch = numel(xsims_Cond_S);

    % append outer most list
    X{condsel} = X_Cond;
end

end
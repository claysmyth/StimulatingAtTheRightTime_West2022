function X_resim = DRL_resim_wUOnly(X, Rorg, m, p)
% resimulate state in the absence of external input
% for debugging purposes

% make a local copy of R first
R = deepCopyStruct(Rorg, 1);
R.IntP.intFx = @spm_fx_compile_120319_wStepEnd;                                   % change simulation function handle

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

%     % form arrays for debugging
%     % TODO: remove
%     for i = 1:numEpoch
%         vec_xsimsComp{i} = [X_cond.xsims.Buffer{i}, X_cond.xsims.SPrev{i}, ...
%             X_cond.xsims.S{i}];
%     end
%     vec_xsims_glComp = X_cond.xsims_gl.S;
%     vec_xsims_SPrevComp = X_cond.xsims.SPrev;
%     vec_xsims_SComp = X_cond.xsims.S;

    % now initialize parallelized for loop
    parfor i = 1:numEpoch
        % unpack meta structure
        R_i = vecR{i};
        m_i = vec_m{i};
        p_i = vec_p{i};

        % form arrays necessary for resimulation
        uFull_i = vecUFull{i};                                              % full intrinsic noise
        xsimsBuff_i = xsimsBuff{i};                                         % initial states for simulation
        tStepEnd_i = size(uFull_i{1}, 1);                                   % end time step for simulation
        tvecRel_i = vec_tvecRel{i};                                         % relative time indices for getting right LFP data

        % re-run simulation
        [xsims_resim_i, xsims_gl_resim_i] = ...
            computeSimDataDRL_wUOnly(R_i, xsimsBuff_i, ...
            m_i, uFull_i, p_i, tStepEnd_i, tvecRel_i);

        % append to outer array
        xsims_resim_Cond{i} = xsims_resim_i{condsel};
        xsims_gl_resim_Cond{i} = xsims_gl_resim_i{condsel};
    end

    % perform quick epoching in resimulated state
    idxSPrev = [zeros(size(X_cond.xsims.tvec_Buffer{1})), ...
        ones(size(X_cond.xsims.tvec_SPrev{1})), ...
        zeros(size(X_cond.xsims.tvec_S{1}))];
    idxS = [zeros(size(X_cond.xsims.tvec_Buffer{1})), ...
        zeros(size(X_cond.xsims.tvec_SPrev{1})), ...
        ones(size(X_cond.xsims.tvec_S{1}))];
    for i = 1:numEpoch
        xsims_resim_Cond_SPrev{i} = xsims_resim_Cond{i}(:, logical(idxSPrev));
        xsims_resim_Cond_S{i} = xsims_resim_Cond{i}(:, logical(idxS));
    end

    % append to output structure
    X_resim_Cond.xsims.SPrev = xsims_resim_Cond_SPrev;
    X_resim_Cond.xsims.S = xsims_resim_Cond_S;
    X_resim_Cond.xsims_gl.S = xsims_gl_resim_Cond;

    X_resim{condsel} = X_resim_Cond;

%     % now test equality for debugging purposes
%     for i = 1:numEpoch
%         try
%             assert(allclose(xsims_resim_Cond{i}, vec_xsimsComp{i}));
%             assert(allclose(xsims_gl_resim_Cond{i}, vec_xsims_glComp{i}));
%         catch
%             t1 = 1;
%         end
%     end
% 
%     for i = 1:numEpoch
%         assert(allclose(xsims_resim_Cond_SPrev{i}, vec_xsims_SPrevComp{i}));
%         assert(allclose(xsims_resim_Cond_S{i}, vec_xsims_SComp{i}));
%     end

end
%% pre-ample to all scripts

clear; close all; clc;

%% Initialize paths and configuration structure

R = project_AddPaths();                                                     % organizes file paths and add dependencies to path
R = setupBasalGangliaModelForDRL(R);                                        % add configuration/settings to R

%% Model-based rollout
% Now load model and run simulation without external stim
% generate the state variable
tic;
[R, m, p] = loadModelForDRL(R);                                             % load saved model fitted on rat data
u = simNoiseIntr(R, m);                                                     % obtain the initial intrinsic noise
[~, ~, ~, xsims, xsims_gl, wflag, R] = ...                                     % generate the state variable s
    computeSimData120319(R,m,u,p,0);
timeSim = toc;

%% Perform epoching

% update parameter set R
tic;
R = updateREpochForDRL(R);

% now perform epoching in state and random noise variable
X = DRL_epochData(xsims, xsims_gl, u, R);
timeEpoch = toc;

%% next test with simulation

DRL_resim_wUOnly(X, R, m, p);

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
[xsimsAlt dum wflag] = spm_fx_compile_120319_wStepEnd(...
    R_i, xBuff, uFull, p_i, m_i, tStepEnd);

% Run Observer function
% Subloop is local optimization of the observer gain
glorg = p_i.obs.LF;
gainlist = R_i.obs.glist;
for gl = 1:length(gainlist)
    % run observation function
    p_i.obs.LF = glorg+gainlist(gl);
    if isfield(R_i.obs,'obsFx')
        [xsims_glAlt{gl}, ~, wflag(1)] = ...
            R_i.obs.obsFx(xsimsAlt,m_i,p_i,R_i, tvecRel_i);
    else
        xsims_glAlt{gl} =xsimsAlt;
    end
end
xsims_glAlt = xsims_glAlt{1};

t1 = 1;

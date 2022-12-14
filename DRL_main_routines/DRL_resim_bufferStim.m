function [X_stim, A, R] = DRL_resim_bufferStim(X, Rorg, morg, porg, varargin)
% In current iteration, varargin is an optional parameter for passing in 
% an external stimulation vector. Size(varargin) = 1.

%% Input parsing

% Handle the optional inputs
p = inputParser;
p.KeepUnmatched = true;

addParameter(p, 'exStimStruct', {});                                        % optional handle for passing in external stim

addParameter(p, 'SScomb', 1, ...                                            % stimulation type
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'stimAmp', 1e3, ...                                         % stimulation amplitude
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'minNPart', 2, ...                                          % min # of random stim parts
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'maxNPart', 8, ...                                          % max # of random stim parts
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'minLenPart', 0.1, ...                                      % min length of each stim part in seconds
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

addParameter(p, 'ampScaleFactor', [0.9, 1.1], ...                           % scaling factor for amplitude in random stim
    @(x) validateattributes(x, {'double'}, {'nonempty'}));

parse(p,varargin{:});

% Handles incorrect inputs
UnmatchedParam = fieldnames(p.Unmatched);
if ~isempty(UnmatchedParam)
    error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
end

% unpacking variable
exStimStruct = p.Results.exStimStruct;
SScomb = p.Results.SScomb;
stimAmp = p.Results.stimAmp;

minNPart = p.Results.minNPart;
maxNPart = p.Results.maxNPart;
minLenPart = p.Results.minLenPart;
ampScaleFactor = p.Results.ampScaleFactor;

% convert from seconds to samples
minLenPart = fix(minLenPart / Rorg.IntP.dt);

%% Load in the model

R = deepCopyStruct(Rorg, 1);

for i = 1 % These are the different stim types
    %% Define stimulation conditions
    if SScomb == 1
        % Stimulating  STN - STN cDBS- DC modulation
        senssite = 4;
        stimsite = 1;
        R.IntP.phaseStim.stimFx = @dcStimDRL_v1;

        % phflag = 0;
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = stimAmp; % 1e4;%
    
    elseif SScomb == 2
        % Stimulating  STN - STN DBS with random DC modulation
        senssite = 4;
        stimsite = 1;
        R.IntP.phaseStim.stimFx = @getRandomStimulationAction;
        % phflag = 0;
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = stimAmp; % 1e4;%
    
    elseif SScomb == 3
        % Stimulating STN - STN DBS with externally (e.g. DRL output) 
        % determined stim vector
        senssite = 4;
        stimsite = 1;
        R.IntP.phaseStim.stimFx = @passExternalStimVector;
        
        R = typeIstimPars_v3(R);
        R.IntP.phaseStim.stimAmp = stimAmp; % 1e4;%
    end
end

R.IntP.phaseStim.sensStm = [senssite stimsite];
R.IntP.phaseStim.eps = Rorg.IntP.phaseStim.eps;

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
    vec_m = deepCopyStruct(morg, numEpoch);
    vec_p = deepCopyStruct(porg, numEpoch);
    
    % unpack external stim struct in the case of external stim
    if isequal(R.IntP.phaseStim.stimFx, @passExternalStimVector)
        A_cond = exStimStruct{condsel};
    end

    % further unpack for parfor
    for i = 1:numEpoch
        % unpack metadata structure
        R_i = vecR{i};
        p_i = vec_p{i};
        m_i = vec_m{i};
    
        % form necessary variables
        vecUFull{i} = {[X_cond.u.Buffer{i}; X_cond.u.SPrev{i}; ...
            X_cond.u.S{i}]};

        % normalize the intrinsic noise first 
        % (note a single epoch)
        us = vecUFull{i}{1};
        for j = 1:m_i.m
            C    = exp(p_i.C(j));
            us(:,j) = C*us(:,j).*0.01;
            % uexs(:,j) = C*uexs(:,i)*0.01;
        end
        uvar = std(us(:,R_i.IntP.phaseStim.sensStm(2)));
        assert(uvar ~= 0)
        
        % make the stim variables
        if R_i.IntP.phaseStim.switch
            % make actual stim pattern
            if isequal(R_i.IntP.phaseStim.stimFx, @dcStimDRL_v1)
                [uexsS_i, R_i] = R_i.IntP.phaseStim.stimFx(R_i, R_i.IntP.phaseStim.stimAmp, ...
                    size(X_cond.u.S{i}), uvar, ampScaleFactor);
            elseif isequal(R_i.IntP.phaseStim.stimFx, @getRandomStimulationAction)
                [uexsS_i, R_i] = R_i.IntP.phaseStim.stimFx(R_i, ...
                    R_i.IntP.phaseStim.stimAmp, size(X_cond.u.S{i}), ...
                    uvar, minNPart, maxNPart, minLenPart, ampScaleFactor);
            elseif isequal(R_i.IntP.phaseStim.stimFx, ...
                    @passExternalStimVector)
                A_i = A_cond.uexs.S{i};
                [uexsS_i, R_i] = R_i.IntP.phaseStim.stimFx(R_i, ...
                    R_i.IntP.phaseStim.stimAmp, size(X_cond.u.S{i}), ...
                    uvar, A_i);
            end

            % make the stim pattern in the burn in period
            [uexsSPrev_i, ~] = dcStimDRL_v1(R_i, R_i.IntP.phaseStim.stimAmp, ...
                size(X_cond.u.SPrev{i}), uvar, [1, 1]);
            
            % start with initial ramping
            bufferLength = fix(minLenPart * 5);
            SPrevLength = size(uexsSPrev_i, 1);
            svEnd = (SPrevLength - bufferLength):SPrevLength;
            ramp = linspace(uexsSPrev_i(end, R.IntP.phaseStim.sensStm(2)), ...
                uexsS_i(1, R.IntP.phaseStim.sensStm(2)), numel(svEnd));
            uexsSPrev_i(svEnd, R.IntP.phaseStim.sensStm(2)) = ramp;
            
            % make stim during buffer and combine
            uexsBuff_i = zeros(size(X_cond.u.Buffer{i}));
            uexsFull_i = {[uexsBuff_i; uexsSPrev_i; uexsS_i]};
        else
            uexsFull_i = zeros(size(us));
        end
        uexsFull_Cond{i} = uexsFull_i;
        uexsS_Cond{i} = uexsS_i / uvar;

    end
    vec_tvecRel = X_cond.xsims_gl.tvecRel;
    xsimsBuff = X_cond.xsims.Buffer;

%     % form arrays for debugging
%     % TODO: remove
%     vec_xsims_glComp = X_cond.xsims_gl.S;
%     vec_xsims_SPrevComp = X_cond.xsims.SPrev;
    
    % now initialize the parallelized for loop
    parfor i = 1:numEpoch
        % unpack meta structure
        R_i = vecR{i};
        m_i = vec_m{i};
        p_i = vec_p{i};

        % form arrays necessary for resimulation
        uFull_i = vecUFull{i};
        xsimsBuff_i = xsimsBuff{i};                                         % initial states for simulation
        tStepEnd_i = size(uFull_i{1}, 1);                                   % end time step for simulation
        tvecRel_i = vec_tvecRel{i};
        uexsFull_i = uexsFull_Cond{i};

        % re-run simulation w/ DBS Stim
        [xsims_stim_i, xsim_gl_stim_i] = ...
            computeSimDataDRL(R_i, xsimsBuff_i, m_i, uFull_i, uexsFull_i, ...
            p_i, tStepEnd_i, tvecRel_i);

        % append to outer array
        xsims_stim_Cond{i} = xsims_stim_i{condsel};
        xsims_gl_stim_Cond{i} = xsim_gl_stim_i{condsel};
    end

%     % perform quick epoching in resimulated state
%     idxSPrev = [zeros(1, size(X_cond.xsims.Buffer{1}, 2)), ...
%         ones(1, size(X_cond.xsims.SPrev{1}, 2)), ...
%         zeros(1, size(X_cond.xsims.S{1}, 2))];
%     idxS = [zeros(1, size(X_cond.xsims.Buffer{1}, 2)), ...
%         zeros(1, size(X_cond.xsims.SPrev{1}, 2)), ...
%         ones(1, size(X_cond.xsims.S{1}, 2))];
% 
%     for i = 1:numEpoch
%         xsims_stim_Cond_SPrev{i} = xsims_stim_Cond{i}(:, logical(idxSPrev));
%         xsims_stim_Cond_S{i} = xsims_stim_Cond{i}(:, logical(idxS));
%     end
% 
%     % append to output structure for next state
%     X_stim_Cond.xsims.SPrev = xsims_stim_Cond_SPrev;
%     X_stim_Cond.xsims.S = xsims_stim_Cond_S;
    X_stim_Cond.xsims_gl.S = xsims_gl_stim_Cond;

    X_stim{condsel} = X_stim_Cond;

    % append to outer structure for action
    A_Cond.uexs.S = uexsS_Cond;
    A_Cond.uexs.Full = uexsS_Cond;

    A{condsel} = A_Cond;

end
end
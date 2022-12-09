function [uexs,R] = dcStimDRL_v1(R, stimLength, xstore, dt, uvar)
% generate the entire stim sequence for DRL simulations 
% input
%   -R: struct storing all hyperparameter
%   -stimAmp: desired stimulation amplitude
%   -stimLength: desired length of stimulation in seconds
%   -dt: time between two samples, 1/fs
%   -uvar: normalization constant?


% Sanity check: if threshold is not set yet then return NaN
if R.IntP.phaseStim.eps == 0
    R.IntP.phaseStim.eps = nan;
    return
end

% calculate the stim amplitude
A = (stimAmp * uvar);

% create empty array for total stimulation
uexs = zeros(size(xstore))


if (tstep+fix(R.IntP.phaseStim.stimlength/dt))<size(uexs,1) || tstep==0

    
    % BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
    if R.IntP.phaseStim.eps == 0
        R.IntP.phaseStim.eps = nan;
        return
    end
    
    %  ts = 0:dt:3;
    A = (R.IntP.phaseStim.stimAmp*uvar); % setup the amplitude of the stim
    off = 0; %uexs(tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) == 0;
    
    if  ~off; % && gate
        sv = tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt);
        uexs(sv,R.IntP.phaseStim.sensStm(2)) = repmat(A,size(sv));
    end
end

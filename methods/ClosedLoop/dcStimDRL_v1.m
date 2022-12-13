function [uexs,R] = dcStimDRL_v1(R, stimAmp, stimSize, uvar)
% generate the entire stim sequence for DRL simulations 
% input
%   -R: struct storing all hyperparameter
%   -stimAmp: desired stimulation amplitude
%   -stimSize: desired shape of stimulation
%   -uvar: normalization constant?


% Sanity check: if threshold is not set yet then return NaN
if R.IntP.phaseStim.eps == 0
    R.IntP.phaseStim.eps = nan;
    return
end

% calculate the stim amplitude
A = (stimAmp * uvar);

% create empty array for total stimulation
uexs = zeros(stimSize);
uexs(:, R.IntP.phaseStim.sensStm(2)) = A;


end

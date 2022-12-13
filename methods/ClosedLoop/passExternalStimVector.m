function [uexs,R] = passExternalStimVector(R, stimAmp, stimSize, uvar, A_i)
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

% also double check the shape of stimulation
assert(all(size(A_i) == stimSize))

% scale the external input
uexs = A_i * uvar;

end
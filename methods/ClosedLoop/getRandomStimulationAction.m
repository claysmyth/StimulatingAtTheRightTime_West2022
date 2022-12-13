function [uexs, R] = getRandomStimulationAction(R, stimAmp, stimSize, uvar, ...
    minNPart, maxNPart, minLenPart, ampScaleFactor)
% generate the entire stim sequence for DRL simulations 
% input
%   -R: struct storing all hyperparameter
%   -stimAmp: desired stimulation amplitude
%   -stimSize: desired shape of stimulation
%   -uvar: normalization constant?
%   -minNPart: min # of partititions
%   -maxNPart: max # of partititions
%   -minLenPart: min length of each partititions in samples
%   -ampScaleFactor: range of scaling of stimulation


% Sanity check: if threshold is not set yet then return NaN
if R.IntP.phaseStim.eps == 0
    R.IntP.phaseStim.eps = nan;
    return
end

num_partitions = round(unifrnd(minNPart, maxNPart));
partition_inds = [0, 0]; % Dummy instantiation

action_length = stimSize(1);

while true
    partition_inds = sort(randi(action_length, 1, num_partitions-1)); % Randomly partition action sequence
    partition_inds = [1, partition_inds, action_length];

    if all(diff(partition_inds) > minLenPart)
        break
    end
end

%  ts = 0:dt:3;
A = (stimAmp * uvar); % setup the amplitude of the stim
uexs = zeros(stimSize);
bufferLength = fix(minLenPart);

for i=1:length(partition_inds)-1
    sv = partition_inds(i):partition_inds(i+1);
    random_scale_factor = unifrnd(ampScaleFactor(1), ampScaleFactor(2)); % For randomly scaling the dc stim amplitude
    currAmp = A * random_scale_factor;
    uexs(sv,R.IntP.phaseStim.sensStm(2)) = currAmp;

    % also make the ramping up
    if i > 1
        svBuffer = (partition_inds(i) - bufferLength):partition_inds(i);
        prevAmp = uexs(partition_inds(i) - 1, R.IntP.phaseStim.sensStm(2));
        ramp = linspace(prevAmp, currAmp, numel(svBuffer));
        uexs(svBuffer, R.IntP.phaseStim.sensStm(2)) = ramp;

    end

end

end

function [uexs,R, stim_action] = getRandomStimulationAction(R, state_length, policy_case, buffer_length, max_number_partitions)
% Returns a random stim_action ( 1 x (state_length+buffer_length) ) pulled according to policy_case
stim_action = zeros(state_length)
partition_inds = [0, 0]

if strcmp(policy_case, 'sinusoid_bursts') | strcmp(policy_case, 'random_bursts')
    num_partitions = unifrnd(2, max_number_partitions)
    while min(diff(partitions_inds)) < min_partition_size
        partition_inds = sort(randi(state_length, 1, num_partitions-1))
    partitions_to_stim = sort(randsample(num_partitions, randi(num_partitions)))

switch policy_case:
    case 'cDBS'
        % stim_action is a single scalar value for the entire vector.
        random_scale_factor = unifrnd(-2,2)
        stim_action = (stim_action + 1) * (R.IntP.phaseStim.stimAmp*random_scale_factor*uvar)
    case 'sinusoid_bursts'
        % stim action is a random assortment of sinusoid bursts. Sinusoid frequencies are chosen guassian around 18 Hz
        while incomplete:
            random_scale_factor = unifrnd(0.1,2)
            random_phase
            imcomplete = false
        end
    case 'random_bursts'
        while incomplete:

% Integrate stim_action into uexs below
sv = tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt);
stim = zeros(size(sv));
stim(zci) = A;
stim(:) = sum(stim(zci))/numel(stim);
uexs(sv,R.IntP.phaseStim.sensStm(2)) = stim;

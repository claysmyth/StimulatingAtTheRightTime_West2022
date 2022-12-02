function [uexs,R,phi] = getRandomStimulationAction(uexs,R,tstep,xstore,dt,uvar,phi, min_number_partitions, max_number_partitions, min_partition_length, amplitude_range)
if (tstep+fix(R.IntP.phaseStim.stimlength/dt))<size(uexs,1) || tstep==0

    
    % BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
    if R.IntP.phaseStim.eps == 0
        R.IntP.phaseStim.eps = nan;
        return
    end

    num_partitions = unifrnd(min_number_partitions, max_number_partitions);
    partition_inds = [0, 0]; % Dummy instantiation

    action_length = tstep+fix(R.IntP.phaseStim.stimlength/dt);

    while min(diff(partitions_inds)) < min_partition_length
        partition_inds = sort(randi(action_length, 1, num_partitions-1)); % Randomly partition action sequence
    end

    partition_inds = partition_inds + tstep; %shift partition inds to align with timeseries
    
    %  ts = 0:dt:3;
    A = (R.IntP.phaseStim.stimAmp*uvar); % setup the amplitude of the stim
    off = 0; %uexs(tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) == 0;
    
    if  ~off; % && gate
        for i=1:length(partition_inds)-1
            sv = partition_inds(i):partition_inds(i+1);
            random_scale_factor = unifrnd(amplitude_range(1), amplitude_range(2)); % For randomly scaling the dc stim amplitude
            uexs(sv,R.IntP.phaseStim.sensStm(2)) = repmat(A * random_scale_factor,size(sv)); % Insert random amplitude DC stim into stimulation struct
            phi(sv,R.IntP.phaseStim.sensStm(2)) = 0;
        end
    end
end
% Returns a random stim_action ( 1 x (state_length+buffer_length) ) pulled according to policy_case
% stim_action = zeros(state_length)
% partition_inds = [0, 0]
% 
% if strcmp(policy_case, 'sinusoid_bursts') | strcmp(policy_case, 'random_bursts')
%     num_partitions = unifrnd(2, max_number_partitions)
%     while min(diff(partitions_inds)) < min_partition_size
%         partition_inds = sort(randi(state_length, 1, num_partitions-1))
%     partitions_to_stim = sort(randsample(num_partitions, randi(num_partitions)))
% 
% switch policy_case:
%     case 'cDBS'
%         % stim_action is a single scalar value for the entire vector.
%         random_scale_factor = unifrnd(-2,2)
%         stim_action = (stim_action + 1) * (R.IntP.phaseStim.stimAmp*random_scale_factor*uvar)
%     case 'sinusoid_bursts'
%         % stim action is a random assortment of sinusoid bursts. Sinusoid frequencies are chosen guassian around 18 Hz
%         while incomplete:
%             random_scale_factor = unifrnd(0.1,2)
%             random_phase
%             imcomplete = false
%         end
%     case 'random_bursts'
%         while incomplete:
% 
% % Integrate stim_action into uexs below
% sv = tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt);
% stim = zeros(size(sv));
% stim(zci) = A;
% stim(:) = sum(stim(zci))/numel(stim);
% uexs(sv,R.IntP.phaseStim.sensStm(2)) = stim;

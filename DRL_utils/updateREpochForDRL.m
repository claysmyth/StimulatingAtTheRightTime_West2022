function R = updateREpochForDRL(R)
% update relevant parameters in R for epoching for DRL

% update observe function to avoid warning with brnIn time
R.obs.obsFx = @DRL_observe_data;

% first determine the number of epochs based on brn time and
% reset sim time
R.DRL.brnNEpoch = ceil(R.obs.brn / R.DRL.epLen);                            % number of epochs per burn in time
R.DRL.brnNSample = R.obs.brn / R.IntP.dt;                                   % number of samples in epochs to discard
R.DRL.tendReCalc = floor(R.IntP.tend / R.DRL.epLen) * R.DRL.epLen; 
if R.DRL.tendReCalc ~= R.IntP.tend                         
    R.IntP.tend = R.DRL.tendReCalc;                                             
    R.IntP.nt = R.IntP.tend/R.IntP.dt;      
    R.IntP.tvec = linspace(0,R.IntP.tend,R.IntP.nt);
    R.IntP.tvec_obs(R.IntP.tvec_obs > R.IntP.tend) = [];
end

end
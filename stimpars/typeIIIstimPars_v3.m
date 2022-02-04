function R = typeIIIstimPars_v3(R)
R.IntP.phaseStim.filtflag = 0;
R.IntP.phaseStim.buff = 1; % This is the buffer used to compute the current phase
R.IntP.phaseStim.minBS = 5/1000; %  Minimum burst length
R.IntP.phaseStim.trackdelay = 25/1000; % this is the delay to take the signal (current tstep-delay)
R.IntP.phaseStim.upperiod  = 25/1000; % update period
R.IntP.phaseStim.stimlength = 0.3; % 300ms stim delivery
R.IntP.phaseStim.stimAmp = 1/4; % times the variance of the normal input;
R.IntP.phaseStim.stimPeriod = 0.5;  % stimulation period
R.IntP.phaseStim.stimGap = 1; % break between stimulation bouts

R.IntP.phaseStim.filtflag = 0;
R.IntP.phaseStim.epsthresh = 75; 
R.IntP.phaseStim.eps = 0;
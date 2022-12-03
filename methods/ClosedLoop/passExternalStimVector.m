function [uexs,R,phi] = passExternalStimVector(uexs,R,tstep,xstore,dt,uvar,phi, external_stim_vector)
if (tstep+fix(R.IntP.phaseStim.stimlength/dt))<size(uexs,1) || tstep==0

    
    % BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
    if R.IntP.phaseStim.eps == 0
        R.IntP.phaseStim.eps = nan;
        return
    end
    
    %  ts = 0:dt:3;
    %A = (R.IntP.phaseStim.stimAmp*uvar); % setup the amplitude of the stim
    off = 0; %uexs(tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) == 0;
    
    if  ~off; % && gate
        %sv = tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt);
        sv = tstep:tstep+size(external_stim_vector);
        uexs(sv,R.IntP.phaseStim.sensStm(2)) = external_stim_vector;
        phi(sv,R.IntP.phaseStim.sensStm(2)) = 0;
    end
end
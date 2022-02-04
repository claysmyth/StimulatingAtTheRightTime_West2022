function [uexs,R,phi] = highFreqStim_pulse_v1(uexs,R,tstep,xstore,dt,uvar,phi)
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
        phiPred = 2*pi*130.*linspace(0,size(sv,2)*dt,size(sv,2));
        zci = find(diff(sign(sin(phiPred)))>1);
%         zci = [zci find(diff(sign(sin(phiPred)))>1)+1];

        stim = zeros(size(sv));
        stim(zci) = A;
        uexs(sv,R.IntP.phaseStim.sensStm(2)) = stim;
        phi(sv,R.IntP.phaseStim.sensStm(2)) = phiPred;
    end
end

demo = 0;
% Demo only
% These are for demo only
if demo
    BA = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(2)),:);
    BUp = padarray(BA,[0 1/dt]);
    BUB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,BUp);
    BUB([1:1/dt 1+end-1/dt:end]) = [];
    
    % Hilbert
    BUBp = padarray(BUB,[0 1/dt]); % pad again
    BAEnv = abs(hilbert(BUBp));
    BPhi = angle(hilbert(BUBp));
    
    % Remove Padding
    BAEnv([1:1/dt end-1/dt:end]) = [];
    BPhi([1:1/dt end-1/dt:end]) = [];
    
    
    figure(1)
    clf
    plot(BEnv(end-fix(0.3/dt)-fix(R.IntP.phaseStim.minBS/dt):end-fix(R.IntP.phaseStim.trackdelay/dt)))
    hold on
    plot([0 numel(uexs)],[R.IntP.phaseStim.eps R.IntP.phaseStim.eps],'k--')
    xlim([0 600])
    drawnow
    a = 1;
    figure(2)
    clf
    yyaxis left; plot(BAEnv); hold on;plot(BA); plot(BUB); yyaxis right; plot(uexs(:,R.IntP.phaseStim.sensStm(2)))
    yyaxis left; plot([0 numel(uexs)],[R.IntP.phaseStim.eps R.IntP.phaseStim.eps],'k--')
    xlim([tstep-2e3 tstep+2e3])
    ylim([-1.5 1.5]*1e-7);
    drawnow
end
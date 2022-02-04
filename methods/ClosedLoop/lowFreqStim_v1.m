function [uexs,R,phi] = lowFreqStim_v1(uexs,R,tstep,xstore,dt,uvar,phi)
if (tstep+fix(R.IntP.phaseStim.stimlength/dt))<size(uexs,1) || tstep==0
    if R.IntP.phaseStim.eps == 0
        BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),:);
        
    else
        BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),end-((R.IntP.phaseStim.buff)/dt):end);
    end
    
    if R.IntP.phaseStim.filtflag == 0
        [dum,B,A] = ft_preproc_bandpassfilter(BU, 1/dt, [14 21],4,'but','twopass');
        R.IntP.phaseStim.filtA = A;
        R.IntP.phaseStim.filtB = B;
        R.IntP.phaseStim.filtflag = 1;
    end
    % Bandpass Filter
    BUp = padarray(BU,[0 1/dt]);
    BUB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,BUp);
    BUB([1:1/dt 1+end-1/dt:end]) = []; 
    
    % Hilbert
    BUBp = padarray(BUB,[0 1/dt]); % pad again
    BEnv = abs(hilbert(BUBp));
    
    % Remove Padding
    BEnv([1:1/dt end-1/dt:end]) = [];
    
    % BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
    if R.IntP.phaseStim.eps == 0
        R.IntP.phaseStim.eps = prctile(BEnv,R.IntP.phaseStim.epsthresh);
        return
    end
    
    %  ts = 0:dt:3;
    A = (R.IntP.phaseStim.stimAmp*uvar); % setup the amplitude of the stim
    % Check the Envelope for gating
    % This means the envelope must be suprathreshold for:
    % t-delay-minBS:t-delay
    gate = all(BEnv(end-fix(R.IntP.phaseStim.trackdelay/dt)-fix(R.IntP.phaseStim.minBS/dt):end-fix(R.IntP.phaseStim.trackdelay/dt)) > R.IntP.phaseStim.eps);
%     instfreq = 18; %median(1./(diff(unwrap(BPhi(end-fix(R.IntP.phaseStim.trackdelay/dt)-fix(R.IntP.phaseStim.minBS/dt)-1:end-fix(R.IntP.phaseStim.trackdelay/dt))))));
    brake = ~all(uexs(tstep-(R.IntP.phaseStim.stimGap/dt):tstep,R.IntP.phaseStim.sensStm(2))==0); % If within
    % This sets up the stim period and ensures breaks
    tv = tstep:tstep+(R.IntP.phaseStim.stimPeriod/dt);
    if gate && ~brake
        uexs(tv,R.IntP.phaseStim.sensStm(2)) = 1e-32;
    end
    off = uexs(tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) == 0;
    
    if  ~off; %  && gate
        sv = tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt);
        phiPred = 2*pi*18.*linspace(0,size(sv,2)*dt,size(sv,2));
        stim = sin(phiPred).*A;
        uexs(sv,R.IntP.phaseStim.sensStm(2)) = stim;
        phi(sv,R.IntP.phaseStim.sensStm(2)) = phiPred;
    end
end
% if brake
%         uexs(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
%     phi(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
% end

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
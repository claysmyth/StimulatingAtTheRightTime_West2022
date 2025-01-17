function [uexs,R,phi] = zeroCrossingPhaseStim_v3(...
    uexs,R,tstep,xstore,dt,uvar,phi,demo)
% This function simulates the phase dependent simulation that simulates at 
% times of zerocrossing of the current phase
%
% output
%   uexs - output external stimulation pattern
%   R - 
%   phi - phase of stimulation?
%
% input
%   uexs - input external stimulation pattern
%   R -  strucutre for experiment variables
%   tstep - ?
%   xstore - simulated local field potential pattern
%   dt - time between two samples
%   uvar - 
%   phi - 
%   demo - flag for whether or not to plot the stimulation pattern
%       simulated

if nargin<8
    demo = 0;
end
if (tstep+fix(R.IntP.phaseStim.stimlength/dt))<size(uexs,1) || tstep==0

    % if eps==0 then extract the entire time series for the sense site
    % otherwise extract a single buffer from the channel
    if R.IntP.phaseStim.eps == 0
        BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),:);
        
    else
        BU = xstore(R.obs.outstates(R.IntP.phaseStim.sensStm(1)),end-((R.IntP.phaseStim.buff)/dt):end);
    end
    
    % designs a bandpass filter in the low beta range
    if R.IntP.phaseStim.filtflag == 0
        [dum,B,A] = ft_preproc_bandpassfilter(BU, 1/dt, [14 21],4,'but','twopass');
        R.IntP.phaseStim.filtA = A;
        R.IntP.phaseStim.filtB = B;
        R.IntP.phaseStim.filtflag = 1;
    end
    
    % If you want to test the effect of observation noise
    if isfield(R.IntP.phaseStim,'obsSNR')
        if ~isinf(R.IntP.phaseStim.obsSNR) % unless SNR is inf (i.e., no obs noise)
            X =(1/(10.^(R.IntP.phaseStim.obsSNR/10)));
            X = X.*std(BU);
            BU = BU + X.*randn(size(BU));
        end
    end
    
    
    
    % Bandpass Filter
    BUp = padarray(BU,[0 1/dt]);
    BUB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,BUp);
    BUB([1:1/dt 1+end-1/dt:end]) = [];
    
    % Hilbert
    BUBp = padarray(BUB,[0 1/dt]); % pad again
    BEnv = abs(hilbert(BUBp));
    BPhi = angle(hilbert(BUBp));
    
    % Remove Padding
    BEnv([1:1/dt end-1/dt:end]) = [];
    BPhi([1:1/dt end-1/dt:end]) = [];
    
    % BEnv = smooth(abs(BUB),200); %abs(HB); %smooth(abs(HB),100);
    % determines the threshold for beta burst events
    if R.IntP.phaseStim.eps == 0
        R.IntP.phaseStim.eps = prctile(BEnv,R.IntP.phaseStim.epsthresh);
        return
    end
    
    %  ts = 0:dt:3;
    A = (R.IntP.phaseStim.stimAmp*uvar); % setup the amplitude of the stim
    % Check the Envelope for gating
    gate = all(BEnv(end-fix(R.IntP.phaseStim.trackdelay/dt)-fix(R.IntP.phaseStim.minBS/dt):end-fix(R.IntP.phaseStim.trackdelay/dt)) > R.IntP.phaseStim.eps);
    instfreq = 18; % Fixed
    %     instfreq = median(1./(diff(unwrap(BPhi(end-fix(R.IntP.phaseStim.trackdelay/dt)-fix(R.IntP.phaseStim.minBS/dt)-1:end-fix(R.IntP.phaseStim.trackdelay/dt))))));
    brake = ~all(uexs(tstep-(R.IntP.phaseStim.stimGap/dt):tstep,R.IntP.phaseStim.sensStm(2))==0); % If within
    % This sets up the stim period and ensures breaks. Stimulation
    % is delivered:
    if gate && ~brake
        uexs(tstep:tstep+(R.IntP.phaseStim.stimPeriod/dt),R.IntP.phaseStim.sensStm(2)) = 1e-32;
    end
    off = uexs(tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) == 0;
    
    if  ~off  %&& gate
        zci = find(diff(sign(BUB))>1,1,'last'); % location of last positive zero crossing
        zstep = (length(BUB)-zci); % location relative to curren time step
        
        % zstep*dt is the current time of the sinusoid relative to t = 0;
        cph = zstep*dt;
        phiPred = 2*pi*instfreq*(cph:dt:cph+R.IntP.phaseStim.stimlength)- pi/2;
        
        stim = sin(phiPred + R.IntP.phaseStim.phaseshift).*A;
        uexs(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = stim;
        phi(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = phiPred;
    end
end
% if brake
%         uexs(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
%     phi(tstep:tstep+fix(R.IntP.phaseStim.stimlength/dt),R.IntP.phaseStim.sensStm(2)) = 0;
% end

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
    yyaxis left;cla; plot(BAEnv); hold on;plot(BA); plot(BUB);
    yyaxis right;cla; plot(uexs(:,R.IntP.phaseStim.sensStm(2)))
    yyaxis left; plot([0 numel(uexs)],[R.IntP.phaseStim.eps R.IntP.phaseStim.eps],'k--')
    xlim([tstep-5e2 tstep+5e2])
    %     ylim([-1.5 1.5]*1e-5);
    drawnow
end
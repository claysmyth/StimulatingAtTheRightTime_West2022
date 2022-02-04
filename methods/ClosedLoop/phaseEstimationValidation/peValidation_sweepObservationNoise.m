clear; close all
% DEMO AND VALIDATIONOF PHASE ESTIMATION ALGORITHM
R.IntP.phaseStim.sensStm = [1 1];
R.IntP.dt = 5e-4;
R.IntP.buffer = 50;
R.obs.brn = 3;
R.obs.outstates = 1;


R.IntP.phaseStim.switch = 1;
R = typeIstimPars_v3(R);
R.IntP.phaseStim.epsthresh = 0;
R.IntP.phaseStim.phaseshift = 0;
R.IntP.phaseStim.stimGap = 0;

obsNVec = linspace(-20,10,7);

%% Start some tests
% TEST (1) = NOISY SINE
% TEST (2) = SIMULATED DATA
% demo = 1/0 animation on or off
demo = 0;
for test = 1
     R.IntP.phaseStim.filtflag = 0;
    if test == 2
        load exampledata
    elseif test == 1
        tData=80.0005;  % in Secs
        SR=1/R.IntP.dt;
        nData=SR*tData + 1;
        tmaxis=(0:nData-1)*R.IntP.dt;
        
        fr1=18;
        am1=1;
        am2=2;
        
        dt1=am1*sin(2*pi*fr1*tmaxis+0);
        dt2=am2*(rand(1,nData)-0.5);
        dta=dt1+dt2;
        tmpdata = dta;
    end
    for obsN = 1:numel(obsNVec)
    %% add the noise
    R.IntP.phaseStim.obsSNR = 10; %obsNVec(obsN);
    

            rl = 25;

        epsStim = zeros(numel(tmpdata),1); dt = R.IntP.dt; uexs = zeros(numel(tmpdata),1); phi = zeros(numel(tmpdata),1);
        for tstep = R.IntP.buffer:numel(tmpdata)
            X = tmpdata(1:tstep);
            
            
            if tstep >((R.obs.brn)/dt) && (rem(tstep,rl) == 0) %&& ~any(uexs(tstep,:))
                if R.IntP.phaseStim.switch
                 [uexs,R,phi] = zeroCrossingPhaseStim_v3(uexs,R,tstep,X,dt,std(tmpdata),phi,demo);
                end
                sprintf('%0.3f',tstep/numel(tmpdata))
            end
        end
        
        % Get true phase
        BU = tmpdata;
        BUB = filtfilt(R.IntP.phaseStim.filtB,R.IntP.phaseStim.filtA,padarray(BU,[0 1/dt]));
        BUB([1:1/dt 1+end-1/dt:end]) = [];
        
        
        ampT = abs(hilbert(BUB));
        highAmpInd = find(ampT>prctile(ampT,50));
        phiT = angle(hilbert(BUB));
        
        pred = wrapToPi(phi(1:numel(phiT)))';
        act = phiT;
        
        PLV(test,obsN) = abs(sum(exp(i*(pred-act))))./numel(act);
        PLV_ha(test,obsN) = abs(sum(exp(i*(pred(highAmpInd)-act(highAmpInd)))))./numel(act(highAmpInd));
        
    end
end


figure
plot(obsNVec(1,:),PLV(1,:),'k')
xlabel('Observation Noise log dB'); ylabel('Accuracy of phase recovery (PLV)')
grid on; box off; axis square

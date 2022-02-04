function [R] = BAA_sim_phaseLockedStim_SICompObsNoise(Rorg)
close all
ip = 0;
% Comopute simulations by sweeping across data
% [R,m,permMod] = getSimModelData_v3(Rorg,modID,simtime);
% mkdir([Rorg.rootn 'data\ModelFit\'])
% save([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')

% OR Load it in:
load([Rorg.rootn 'data\modelfit\SimModelData_M10.mat'],'R','m','p')
par = p;
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
warning('Loading Preloaded model, cant change simtime or model choice!!!')
SScomb = 1;
obsNVec = linspace(-20,10,7);
for obsN = 1:numel(obsNVec)
    
    if  SScomb == 1
        % Stimualating M2
        senssite = 4; % STN
        stimsite = 1; % M2
        stimFx = @zeroCrossingPhaseStim_v3;
        stim_sens = 'stimM2_sensSTN';
        
    end
    phflag = 1;
    % Setup stim parameterss
    R = typeIstimPars_v3(R);
    R.IntP.phaseStim.sensStm = [senssite stimsite];
    R.IntP.phaseStim.stimGap = 0.2;
    R.IntP.phaseStim.stimperiod =1;
    %             R.IntP.phaseStim.stimAmp = 0.5;
    % Simulation Coniditions
    R.obs.csd.df = 0.5;
    R = setSimTime(R,48);
    
    % Trans Options
    %     R.obs.SimOrd = 10;
    R.obs.trans.norm = 0;
    R.obs.gainmeth = {};
    
    % Give all timeseries the same input - makes comparable
    rng(5453)
    m.uset.p.scale = m.uset.p.scale;
    uc = innovate_timeseries(R,m);
    uc{1} = uc{1}.*sqrt(R.IntP.dt);
    XBase =par;
    
    % Phase To be Tested
    R.IntP.intFx = @spm_fx_compile_120319_stim;
    
    if phflag
        phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
        phaseShift = phaseShift(1:12); %12
    else
        phaseShift = 0;
    end
    
    %% Loop through Connections
    CON = 1;
    state = 1
    %% Setup Base Model
    Pbase = XBase;
    
    % Simulate Base Model
    uc_ip{1} = uc;
    R.IntP.phaseStim.switch = 0 ;
    R.IntP.phaseStim.phaseshift = 0;
    R.frqz = 2:0.2:150;
    R.IntP.compFx = @nullComp;
    [~,~,feat_sim_base{1},xsim_gl,xsim_ip_base{1},wflag] = computeSimData120319(R,m,uc_ip{1},Pbase,0);
    if wflag == 1
        a = 1;
    end
    % Work out the threshold
    R.IntP.phaseStim.eps = 0;
    [~,R] = zeroCrossingPhaseStim_v3([],R,0,xsim_gl{1},R.IntP.dt);
    
    % work out base beta
    baseFeat = squeeze(feat_sim_base{1}(1,4,4,1,:));
    bp_base(1) = bpSpec(R.frqz,baseFeat,[14 21]);
    bp_base(2) = bpSpec(R.frqz,baseFeat,[21 30]);
    
    
    %% Do Stimulation
    % modulate the observation SNR
    R.IntP.phaseStim.obsSNR = obsNVec(obsN);
    
    R.IntP.phaseStim.switch = 1;
    m = m; % initialise for parfor
    feat_sim_stim = cell(1,12); pU = cell(1,12);
    %             parfor p = 1:numel(phaseShift)
    parfor p = 1:numel(phaseShift)
        Rpar = R;
        % Modulate the phase
        Rpar.IntP.phaseStim.phaseshift = phaseShift(p);
        Rpar.IntP.phaseStim.stimFx = stimFx;
        % Simulate with Stimulation
        [~,~,feat_sim_stim{p}]  = computeSimData120319(Rpar,m,uc_ip{1},Pbase,0);
        uexs = load([Rpar.rootn 'data\phaseStimSave\stim_tmp_' sprintf('%3.f',1000*Rpar.IntP.phaseStim.phaseshift)],'uexs');
        pU{p} = uexs.uexs(Rpar.IntP.phaseStim.sensStm(2),round(Rpar.obs.brn*(1/R.IntP.dt))+1:end);
        disp([obsN p])
    end
    
    for p = 1:numel(phaseShift)
        phaseFeat(:,p) = squeeze(feat_sim_stim{p}(1,4,4,1,:));
        bp = bpSpec(R.frqz,phaseFeat(:,p),[14 21]);
        deltaBP(1,p) = 100*(bp-bp_base(1))./bp_base(1);
        
        bp = bpSpec(R.frqz,phaseFeat(:,p),[21 30]);
        deltaBP(2,p) = 100*(bp-bp_base(2))./bp_base(2);
    end
    deltaBPStore{obsN} = deltaBP;
    phaseFeatStore{obsN} = phaseFeat;
    
    ip = ip+1;
end


%% Now plot
% Plot Spectra
% subplot(1,3,1)
% plot(R,phaseFeatStore{1}
ARCStat = [];
for obsN = 1:numel(obsNVec)
    % Plot the ARCs
%     if obs 
%     plot(phaseShift,deltaBPStore{obsN})
    hold on
    ARCStat(:,1,obsN) = [min(deltaBPStore{obsN}(1,:)) max(deltaBPStore{obsN}(1,:))];
    ARCStat(:,2,obsN) = [min(deltaBPStore{obsN}(2,:)) max(deltaBPStore{obsN}(2,:))];
end
subplot(1,2,1)
plot(obsNVec,squeeze(ARCStat(2,1,:)),'r')
hold on
plot(obsNVec,squeeze(ARCStat(1,1,:)),'b')
xlabel('sensing SNR (log dB)'); ylabel('beta modulation (%)')
grid on; box off; axis square
legend({'Max. ARC','Min. ARC'})

subplot(1,2,2)
plot(obsNVec,squeeze(ARCStat(2,2,:)),'r')
hold on
plot(obsNVec,squeeze(ARCStat(1,2,:)),'b')
xlabel('sensing SNR (log dB)'); ylabel('beta modulation (%)')
grid on; box off; axis square
legend({'Max. ARC','Min. ARC'})
a = 1;
% ! shutdown /s
function bp = bpSpec(Hz,Pxy,flim)
bandInds = find(Hz>=flim(1) & Hz<= flim(2));
bp = sum(Pxy(bandInds))*(numel(bandInds)*(Hz(2)-Hz(1)));
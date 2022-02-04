function [R] = BAA_sim_SI_StimSweep(Rorg)
close all
ip = 0;
% Comopute simulations by sweeping across data
% [R,m,permMod] = getSimModelData_v3(Rorg,modID,simtime);
% mkdir([Rorg.rootn 'data\ModelFit\'])
% save([Rorg.rootn 'data\ModelFit\SimModelData.mat'],'R','m','permMod')

% OR Load it in:
load([Rorg.rootn 'data\modelfit\SimModelData_M10.mat'],'R','m','permMod')
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
warning('Loading Preloaded model, cant change simtime or model choice!!!')

threshList = linspace(0,95,16);
ampList = linspace(-20,20,16);
perList = linspace(0.1,5,16);
for parI = 1:16
    
    %% Define stimulation conditions
    
    % Stimulating  STN
    senssite = 4; % M2
    stimsite = 4; % STN
%             stimFx = @highFreqStim_integrated_v1;
    stimFx = @highFreqStim_pulse_v1;
    phflag = 0;
    R = typeIIIstimPars_v3(R);
    %         R.IntP.phaseStim.stimAmp = ampList(parI);
    %         R.IntP.phaseStim.epsthresh = threshList(parI);
    R.IntP.phaseStim.stimGap = 0;
    R.IntP.phaseStim.stimPeriod = 1;
    R.IntP.phaseStim.sensStm = [senssite stimsite];
    
    
    %     % Simulation Coniditions
    R.obs.csd.df = 0.5;
    R = setSimTime(R,32);
    
    % Trans Options
    %     R.obs.SimOrd = 10;
    R.obs.trans.norm = 0;
    R.obs.gainmeth = {};
    
    % Observe Middle layers
    R.obs.outstates(1) = 3; % change to middle layer
    m.outstates{1} = [0 0 1 0 0 0 0 0];
    
    % Give all timeseries the same input - makes comparable
    rng(5453)
    m.uset.p.scale = m.uset.p.scale;
    uc = innovate_timeseries(R,m);
    uc{1} = uc{1}.*sqrt(R.IntP.dt);
    XBase = permMod{1}.par_rep{1};
    
    % Phase To be Tested
    R.IntP.intFx = @spm_fx_compile_120319_stim;
    
    
    %% Setup Base Model
    Pbase = XBase;
    
    % Simulate Base Model
    uc_ip{1} = uc;
    R.IntP.phaseStim.switch = 0 ;
    R.IntP.phaseStim.phaseshift = 0;
    R.frqz = 2:0.2:150;
    R.IntP.compFx = @nullComp;
    [~,~,feat_sim_base{1},xsim_gl,xsim_ip_base{1},~,Rout] = computeSimData(R,m,uc_ip{1},XBase,0);
    % Work out the threshold
    
    % Compute base power
    baseFeat = squeeze(feat_sim_base{1}(1,4,4,1,:));
    bp_base(1) = bpSpec(Rout.frqz,baseFeat,[14 21]);
    bp_base(2) = bpSpec(Rout.frqz,baseFeat,[21 30]);
    
    %% Do Stimulation
    R.IntP.phaseStim.switch = 1;
    m = m; % initialise for parfor
    deltaBP = zeros(2,16);
    parfor parJ = 1:16
        Rpar = R;
        Rpar.IntP.phaseStim.stimAmp = ampList(parJ);
        Rpar.IntP.phaseStim.epsthresh = 50; %threshList(parJ);
        Rpar.IntP.phaseStim.stimPeriod = perList(parI);

        % Recompute threshold
        Rpar.IntP.phaseStim.eps = 0;
        [~,Rpar] = zeroCrossingPhaseStim_v3([],Rpar,0,xsim_gl{1},Rpar.IntP.dt);
        
        % Modulate the phase
        Rpar.IntP.phaseStim.phaseshift = 0;
        Rpar.IntP.phaseStim.stimFx = stimFx;
        % Simulate with Stimulation
        [~,~,feat_sim_stim]  = computeSimData(Rpar,m,uc_ip{1},Pbase,0);
        
        phaseFeat = squeeze(feat_sim_stim(1,4,4,1,:));
        bp1 = bpSpec(Rout.frqz,phaseFeat,[14 21]);
        bp2 = bpSpec(Rout.frqz,phaseFeat,[21 30]);
        deltaBP(:,parJ) = [100*(bp1-bp_base(1))./bp_base(1) 100*(bp2-bp_base(2))./bp_base(2)];
        
        disp([parJ parI])
    end
    deltaBPMat(:,:,parI) = deltaBP;
    rmdir([R.rootn 'data\phaseStimSave\'],'s')
    % Create save version
end

% Make figure
frqtit = {'B1','B2'};
for sp = 1:2
    subplot(1,2,sp)
    imagesc(perList,ampList,squeeze(deltaBPMat(sp,:,:)))
    set(gca,'YDir','normal')
    ylabel('Amplitude')
    xlabel('Gating threshold (nth percentile)')
    title(frqtit{sp})
%     caxis([-2 20])
    C = colorbar;
    C.Label.String = '% power modulation';
end

rootan = [Rorg.rootn 'data\stimSweep'];
mkdir(rootan)
save([rootan '\BB_' Rorg.out.tag '_stimSweep_HFInt_StimGap.mat'],'deltaBPMat','R','threshList','ampList')


% ! shutdown /s
function bp = bpSpec(Hz,Pxy,flim)
bandInds = find(Hz>=flim(1) & Hz<= flim(2));
bp = sum(Pxy(bandInds))*(numel(bandInds)*(Hz(2)-Hz(1)));
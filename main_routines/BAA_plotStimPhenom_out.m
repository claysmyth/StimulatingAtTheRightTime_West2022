function [R] = BAA_plotStimTypePhenom_out(Rorg)
% close all

% Setup simulations to retrieve
state = 1; % baseline state only
CON = 2; % dont need to discriminate

% STIM Types:
% [1 - baseline shift model of DBS]

phaseN = [1]; % Number of phases for each stim type

% Counters for plots
SSc = 0;
splist = [1 4 2 5 3 6];


SScomb = 8
%% Loop through Connections
% Load Data
rootan = [Rorg.rootn 'data\phaseLockedStim'];
%     load([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],'feat_sim_save')
load([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_Rout' num2str(SScomb) '.mat'],'R')
load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_xsim' num2str(SScomb) '.mat'],'xsim_ip');
load([rootan '\BB_' Rorg.out.tag '_phaseLockedStim_CON_' num2str(CON) '_pU_save' num2str(SScomb) '.mat'],'pU_save')

%
fsamp = 1/R.IntP.dt;
winsize = [0.075*fsamp; 0.075*fsamp]


% Get bandpower of unstimulated (base) model
Rout = R;
deltaBP = []; phaseFeat = [];
p = 1;
% Get Burst Inds
pU = pU_save{state}{p};
stimSelInds = SplitVec(find(abs(pU)>0),'consecutive'); % Split up data based upon stimulation gating
subplot(1,2,2)
[F,Hz] = pwelch(xsim_ip{1}{1}{1}(4,:)',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
baseFeat = F; % this is the STN power spectra
if p == 1
    po(1) = plot(Hz,baseFeat,'k','LineWidth',1.2); hold on
end
bp_base = bpSpec(Hz,baseFeat,[14 21]);

% Estimate change in B1 bandpower per phase

[F,Hz] = pwelch(xsim_ip{2}{p}{1}(4,:)',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
phaseFeat(:,p) = F; % this is the STN power spectra
bp = bpSpec(Hz,phaseFeat(:,p),[14 21]);
deltaBP(p) = 100*((bp-bp_base)./bp_base); % change in BP from baseline

if numel(deltaBP)>1
    [~,supp] = min(deltaBP);
    [~,amp] = max(deltaBP);
    po(2) = plot(Hz,phaseFeat(:,supp),'b','LineWidth',1.2);
    hold on
    po(3) = plot(Hz,phaseFeat(:,amp),'r','LineWidth',1.2);
    legend(po,{'0' num2str(deltaBP(supp),2) num2str(deltaBP(amp),2)},'box','off')
else
    po(2) = plot(Hz,phaseFeat,'g','LineWidth',1.2);
    legend(po,{'0' num2str(deltaBP(1),2) },'box','off')
end
xlabel('Hz'); ylabel('Power'); xlim([2 48])
axis square; box off; grid on;


function bp = bpSpec(Hz,Pxy,flim)
% this finds the bandpower
bandInds = find(Hz>=flim(1) & Hz<= flim(2));
bp = sum(Pxy(bandInds))*(numel(bandInds)*(Hz(2)-Hz(1)));
function BAA_fingerprintCompare_emp(R)
close all
scmap = [ brewermap(4,'Set1'); brewermap(2,'Set2')];
%% Spontaneous Matrices
% Load Canonical States (predefined)
rootan = [R.rootn 'data\ConnectionSweep'];
load([rootan '\BB_' R.out.tag '_EmpiricalStates_ConnectivityMatrix.mat'],'specMatEmp','connectMatEmp')

connectMat = connectMatEmp;
specMat = specMatEmp;
clear connectMatEmp specMatEmp

% Setup State
statelistSpont = [          % state CON
    1 1;    % 1 Control data
    2 1;    % 2 Lesion Data
    ];

aniNum = [7 8];

% connectMat{band,state,CON}
stateDefPLV = []; stateDefSpec = [];
for state = 1:size(statelistSpont,1)
    for aniN = 1:aniNum(state)
        stateDefPLV(:,:,state,aniN) = [abs(connectMat{1,statelistSpont(state,1),statelistSpont(state,2),aniN}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2),aniN}),...
            angle(connectMat{1,statelistSpont(state,1),statelistSpont(state,2),aniN}+connectMat{2,statelistSpont(state,1),statelistSpont(state,2),aniN})];
        stateDefSpec(:,:,state,aniN) = specMat{1,statelistSpont(state,1),statelistSpont(state,2),aniN};
    end
end
SScomb = 1; stm = 2;
matPLVStim = []; matSPECStim = [];

%% Now Load in Stim Data
rootan = [R.rootn 'data\phaseLockedStim'];
load([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysisSweepState_' num2str(1) '_CON_' num2str(1) '_feat' num2str(SScomb) '.mat'],...
    'connectMat','specMat');
%% Stim Matrices
philistStim = 1:12;
%     cmatStim = []; smatStim = [];
for phi = 1:size(philistStim,2)
    matPLVStim(:,:,phi,1) = [abs(connectMat{1,phi,stm}+connectMat{2,phi,stm}) angle(connectMat{1,phi,stm}+connectMat{2,phi,stm})];
    matSPECStim(:,:,phi,1) = specMat{1,phi,stm};
end

%% Do comparison
nmlist = []; plvR2 = nan(size(stateDefPLV,3),size(philistStim,2),max(aniNum)); specR2 = nan(size(stateDefPLV,3),size(philistStim,2),max(aniNum));
for stateSpont = 1:size(stateDefPLV,3)
    for phiStim = 1:size(philistStim,2)
        A{2} = squeeze(matPLVStim(:,:,phiStim,1));
        A{1} =  normalize(squeeze(matSPECStim(:,:,phiStim,1)));
        for aniN = 1:aniNum(stateSpont)
            
            B{2} = squeeze(stateDefPLV(:,:,stateSpont,aniN));
            B{1} =  normalize(squeeze(stateDefSpec(:,:,stateSpont,aniN)));
            
            % Seperate linear part
            Lpart{1} = A{2}(2:6,2:6); 
            Lpart{2} = B{2}(2:6,2:6); 
            [Lpart{1},Lpart{2}] = remnan(Lpart{1}(:),Lpart{2}(:));
            Lpart{1}(Lpart{1}==2) = [];
            Lpart{2}(Lpart{2}==2) = [];
            Lpart{1}(Lpart{1}==0) = [];
            Lpart{2}(Lpart{2}==0) = [];
            [AC{1},BC{1}] = remnan(A{1}(:),B{1}(:));
            
            Lpart{1} = Lpart{1}-mean(Lpart{1});
            Lpart{2} = Lpart{2}-mean(Lpart{2});
            
            %  RSquare
            plvR2(stateSpont,phiStim,aniN) = rsquare(Lpart{1}(:),Lpart{2}(:));
            specR2(stateSpont,phiStim,aniN) =  rsquare(AC{1}(:),BC{1}(:)); subplot(3,1,3);
        end
        
    end
end

% Now realign based on maximum
for aniN = 1:aniNum(stateSpont)
    X = squeeze(specR2(2,:,aniN)); % Realign by lesion
    Y = squeeze(specR2(1,:,aniN));
    [specR2(2,:,aniN),specR2(1,:,aniN),indshift] = realignMax(X,Y);
    
    X = squeeze(plvR2(2,:,aniN)); % Realign by lesion
    Y = squeeze(plvR2(1,:,aniN));
    [plvR2(2,:,aniN),plvR2(1,:,aniN)] = realignMax(X,Y,indshift);
    
end

figure
plotFingerPrintMatch(squeeze(nanmean(plvR2,3)),squeeze(nanmean(specR2,3)),...
    squeeze(nanstd(plvR2,[],3))./sqrt(aniNum'),squeeze(nanstd(specR2,[],3))./sqrt(aniNum'),...
    scmap(5:6,:),1:2,{'Cont.','Lesion'},0,[2 2],1:2);

plotFingerPrintMatch(squeeze(nanmean(plvR2,3)),squeeze(nanmean(specR2,3)),...
    squeeze(nanstd(plvR2,[],3))./sqrt(aniNum'),squeeze(nanstd(specR2,[],3))./sqrt(aniNum'),...
    scmap(5:6,:),1:2,{'Cont.','Lesion'},1,[2 2],3:4);

function LEG = plotFingerPrintMatch(nmlist,slist,nmS,slS,scmap,netstates,legnames,percflag,splotDim,splotInd)
for L = 1:2
    if L == 2
        lm = nmlist;
        ls = nmS;
    elseif  L == 1
        lm = slist;
        ls = slS;
    end
    
    if percflag ==1
        lm = 100*(lm - mean(lm,2));
        ytit = {'Change in '; 'explained variance'};
        %         ytit = {'Similarity to '; 'Spontaneous (R2)'};
    else
        ytit = {'Similarity to '; 'Spontaneous (R2)'};
    end
    
    
    subplot(splotDim(1),splotDim(2),splotInd(L))
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    phaseShiftDeg = rad2deg(phaseShift);
    p = plot(phaseShiftDeg,lm(netstates,:),'LineWidth',2);
    hold on
    ps = plot(phaseShiftDeg,lm(netstates,:)-ls,'LineWidth',1,'LineStyle','--');
    pu = plot(phaseShiftDeg,lm(netstates,:)+ls,'LineWidth',1,'LineStyle','--');
    
    for pip = 1:numel(p)
        p(pip).Color = scmap(netstates(pip),:);
        ps(pip).Color = scmap(netstates(pip),:);
        pu(pip).Color = scmap(netstates(pip),:);
    end
    xlabel('Stimulation Phase (radians)');
    xlim([0 phaseShiftDeg(end)])
    a = gca;
    a.XTick = 0:90:270;
    grid on; box off; axis square
end
subplot(splotDim(1),splotDim(2),splotInd(1))
title('Mapping to State Spectra')
ylabel(ytit)

subplot(splotDim(1),splotDim(2),splotInd(2))
title('Mapping to State PLV')
ylabel(ytit)

set(gcf,'Position',[   655   239   762   523])

LEG = legend(legnames);
LEG.Orientation = 'horizontal';
LEG.Position = [0.304 0.0634 0.437 0.0568];
LEG.Box = 'off';
LEG.Color = 'none';

function [Xshift,Yshift,ind] = realignMax(X,Y,ind)
if nargin<3
    [~,ind] = max(X);
end
Xshift = circshift(X,floor(numel(X)/2)-ind);
Yshift = circshift(Y,floor(numel(X)/2)-ind);


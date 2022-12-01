function BAA_fingerprintCompare_Sweep(R)
close all
scmap = [ brewermap(4,'Set1'); brewermap(2,'Set2')];
SScomb = 1;
rootan = [R.rootn '/data/ConnectionSweep'];

%% Get Connection Strengths for each CON
load([rootan '/BB_' R.out.tag '_ConnectionSweep_CON_1_ck_1.mat'],'ck_1'); % load connection bands (CON 1 and 2 have same data)

% Setup Spontaneous States
statedef =[1 2;
    1 31;
    2 2;
    2 31];
for stateSpont = 1:4
    % Now Load in  Data
    rootan = [R.rootn '/data/phaseLockedStim'];
    load([rootan '/BB_' R.out.tag '_phaseLockedStim_burstAnalysisSweepState_' num2str(statedef(stateSpont,2)) '_CON_' num2str(statedef(stateSpont,1)) '_feat' num2str(SScomb) '.mat'],...
        'connectMat','specMat');
    
    %No Stim matrices
    phi = 1; stm = 1;
    matPLVNOStim(:,:,stateSpont) = [abs(connectMat{1,phi,stm}+connectMat{2,phi,stm}) angle(connectMat{1,phi,stm}+connectMat{2,phi,stm})];
    matSPECNoStim(:,:,stateSpont) = specMat{1,phi,stm};
end

complist = 1:4;

%% Stim data
for CON = 1:2
    matPLVStim = []; matSPECStim = [];
    for state = 1:31
        % Now Load in Stim Data
        rootan = [R.rootn 'data/phaseLockedStim'];
        load([rootan '/BB_' R.out.tag '_phaseLockedStim_burstAnalysisSweepState_' num2str(state) '_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],...
            'connectMat','specMat');
        
        %% Stim Matrices
        % connectMat{band,phi,stm}
        philistStim = 1:12;
        stm = 2;
        %     cmatStim = []; smatStim = [];
        for phi = 1:size(philistStim,2)
            matPLVStim(:,:,phi,state) = [abs(connectMat{1,phi,stm}+connectMat{2,phi,stm}) angle(connectMat{1,phi,stm}+connectMat{2,phi,stm})];
            matSPECStim(:,:,phi,state) = specMat{1,phi,stm};
        end
    end
    %% Do comparison
    plvR2 = []; specR2 = [];
    for state = 1:31
        for phiStim = 1:size(philistStim,2)
            A{2} = squeeze(matPLVStim(:,:,phiStim,state));
            A{1} =  normalize(squeeze(matSPECStim(:,:,phiStim,state)));
            for stateSpont = 1:size(statedef,1)
                B{2} = squeeze(matPLVNOStim(:,:,stateSpont));
                B{1} =  normalize(squeeze(matSPECNoStim(:,:,stateSpont)));
                
                % Seperate linear and circular parts
                Lpart{1} = A{2}(1:6,1:6); Cpart{1} = A{2}(1:6,7:12); % seperate out linear vs circular parts of A
                Lpart{2} = B{2}(1:6,1:6); Cpart{2} = B{2}(1:6,7:12); % seperate out linear vs circular parts of B
                [Cpart{1},Cpart{2}] = remnan(Cpart{1}(:),Cpart{2}(:));
                [Lpart{1},Lpart{2}] = remnan(Lpart{1}(:),Lpart{2}(:));
                Lpart{1}(Lpart{1}==0) = [];
                Lpart{2}(Lpart{2}==0) = [];
                [AC{1},BC{1}] = remnan(A{1}(:),B{1}(:));
                Lpart{1}(Lpart{1}==2) = [];
                Lpart{2}(Lpart{2}==2) = [];
                Cpart{1}(Cpart{1}==0) = [];
                Cpart{2}(Cpart{2}==0) = [];
                
                plot( AC{1}(:))
                hold on
                plot(BC{1}(:))
                
                % RSqr
                plvR2(stateSpont) = rsquare(Lpart{1}(:),Lpart{2}(:));
                specR2(stateSpont) =  rsquare(AC{1}(:),BC{1}(:));
                combR2(stateSpont) =  mean([ specR2(stateSpont)  plvR2(stateSpont)]);
            end
            plvStateList(state,phiStim,CON) = complist(plvR2(complist)==max(plvR2(complist)));
            specStateList(state,phiStim,CON) = complist(specR2(complist)==max(specR2(complist)));
            combStateList(state,phiStim,CON) = complist(combR2(complist)==max(combR2(complist)));
            
            plvR2Store(:,state,phiStim,CON) = (plvR2(complist));
            specR2Store(:,state,phiStim,CON) = (specR2(complist));
            combR2Store(:,state,phiStim,CON) = (combR2(complist));
        end
    end
end


figure(100)
plotFingerPrintMatch(squeeze(plvR2Store(:,1,:,1)),squeeze(specR2Store(:,1,:,1)),...
    squeeze(combR2Store(:,1,:,1)),scmap,1:4,{'HD-Down','HD-Up','PS-Down','PS-Up'},0,[3 3])
subplot(3,3,1); ylim([-1 1]); subplot(3,3,2); ylim([-1 1]); subplot(3,3,3); ylim([-1 1]);


for i = 1:3
    switch i
        case 1
            AB = specStateList;
        case 2 %
            AB = plvStateList;
        case 3
            AB = combStateList;
            AB_r2 = combR2Store;
    end
    for CON = 1:2
        subplot(3,3,(3*CON) + i)
        [lr,zeroind] = min(abs(ck_1(CON,2:31)-1));
        if (ck_1(CON,zeroind+1)-1)<0
            list = [2:zeroind 1 zeroind+2:31];
        else
            list = [2:zeroind-1 1 zeroind+1:31];
        end
        
        if any(diff(ck_1(CON,list))<0)
            error('List is not continuos!')
        end
        D = squeeze(AB(list,1:end,CON));
        D(zeroind,:) = AB(1,:,CON);
        plotMats(ck_1(CON,list),D,scmap(1:4,:),zeroind)
    end
    
end
set(gcf,'Position',[323          95        1435         883])

a = 1;

function LEG = plotFingerPrintMatch(nmlist,slist,clist,scmap,netstates,legnames,percflag,splotDim)
for L = 1:3
    if L == 2
        lm = nmlist;
    elseif  L == 1
        lm = slist;
        %     elseif L == 3
        %         lm = olist;
    elseif L == 3
        lm = clist;
    end
    
    if percflag ==1
        lm = 100*(lm - mean(lm,2));
        ytit = {'Change in '; 'explained variance'};
    else
        ytit = {'Similarity to '; 'Spontaneous (R2)'};
    end
    
    subplot(splotDim(1),splotDim(2),L)
    phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
    phaseShift = phaseShift(1:12); %12
    phaseShiftDeg = rad2deg(phaseShift);
    p = plot(phaseShiftDeg,lm(netstates,:),'LineWidth',2);
    for pip = 1:numel(p)
        p(pip).Color = scmap(netstates(pip),:);
    end
    xlabel('Stimulation Phase (radians)');
    xlim([0 phaseShiftDeg(end)])
    a = gca;
    a.XTick = 0:90:270;
    grid on; box off; axis square
end
subplot(splotDim(1),splotDim(2),1)
title('Mapping to State Spectra')
ylabel(ytit)

subplot(splotDim(1),splotDim(2),2)
title('Mapping to State PLV')
ylabel(ytit)

subplot(splotDim(1),splotDim(2),3)
title('Overall Mapping')
ylabel(ytit)
set(gcf,'Position',[326 239 1091 523])

LEG = legend(legnames);
LEG.Orientation = 'horizontal';
LEG.Position = [0.304 0.0634 0.437 0.0568];
LEG.Box = 'off';
LEG.Color = 'none';

function a = plotMats(y,AB,cmap,zeroind)
phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12
phaseShift = rad2deg(phaseShift);
imagesc(phaseShift,1:numel(y),AB,'AlphaData',~isnan(AB))
colormap(gca,cmap)
a = gca;
a.XTick = [0 90 180 270];
% a.YTick = y
a.YTick = (1:4:numel(y));
a.YTickLabel = num2cell(round(y(1:4:end)*100,0));
set(a, 'ydir', 'normal');
hold on
plot(phaseShift,repmat(zeroind,size(phaseShift)),'k--')
caxis([1 4]); %xlim([0.5 12.5])
axis square;
grid on




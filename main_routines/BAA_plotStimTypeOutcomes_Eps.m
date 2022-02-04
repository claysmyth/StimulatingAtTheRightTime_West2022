function [R] = BAA_plotStimTypeOutcomes_Eps(Rorg)
state = 1; % baseline state only
CON = 1; % dont need to discriminate
close all
SStype = {'stimM2_sensSTN','stimSTN_sensM2','stimSTN_sensGPe','stimM2_sensSTN','stimSTN_sensM2'};
phaseN = [12 12];
% EpsVec = [0:10:100];
% AmpVec = [0:0.25:2];

cmap{1} = [linspace(1,0,7); zeros(1,7); zeros(1,7)]';
cmap{2} = [zeros(1,7); zeros(1,7); linspace(1,0,7)]';
figure

%% Loop through Connections
for SScomb = 1:2
    Rout = Rorg;
    Rout.frqz = 2:0.2:150;
    rootan = [Rorg.rootn 'data\phaseLockedStimEpsComp'];
    load([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_feat' num2str(SScomb) '.mat'],'feat_sim_save','EpsVec','AmpVec')
    
    suppSave = []; ampSave = [];
    for EpsN = 1:numel(EpsVec)
        for AmpN = 1:numel(AmpVec)
            baseFeat = squeeze(feat_sim_save{1,AmpN,EpsN}{1}(1,4,4,1,:));
            bp_base(1) = bpSpec(Rout.frqz,baseFeat,[14 21]);
            bp_base(2) = bpSpec(Rout.frqz,baseFeat,[21 30]);
            state = 1; %:size(ck_1,2)
            deltaBP = []; phaseFeat = [];
            for p = 1:phaseN(SScomb)
                phaseFeat(:,p) = squeeze(feat_sim_save{2,AmpN,EpsN}{p}(1,4,4,1,:));
                bp = bpSpec(Rout.frqz,phaseFeat(:,p),[14 21]);
                deltaBP(1,p) = 100*(bp-bp_base(1))./bp_base(1);
                bp = bpSpec(Rout.frqz,phaseFeat(:,p),[21 30]);
                
                deltaBP(2,p) = 100*(bp-bp_base(2))./bp_base(2);
                
            end
            % low beta
            [suppSave(1,EpsN,AmpN)] = min(deltaBP(1,:));
            [ampSave(1,EpsN,AmpN)] = max(deltaBP(1,:));
            % high beta
            [suppSave(2,EpsN,AmpN)] = min(deltaBP(2,:));
            [ampSave(2,EpsN,AmpN)] = max(deltaBP(2,:));    end
    end
    %     try
    %         load('shiftedCMAP');
    %     catch
    %         disp('Map not available')
    %     end
    bandName = {'B1','B2'};
    for band = 1:2
        subplot(2,3,band+((SScomb-1)*3))
        plot(EpsVec,squeeze(ampSave(band,:,2)),'r','LineWidth',2)
        hold on
        plot(EpsVec,squeeze(suppSave(band,:,2)),'b','LineWidth',2)
        scatter(70,squeeze(suppSave(band,8,2)),100,'kv','filled','LineWidth',2)
        scatter(70,squeeze(ampSave(band,8,2)),100,'ko','filled','LineWidth',2)
        
        axis square; grid on; box off;
        xlabel([bandName{band} ' Threshold (percentile)']); ylabel([bandName{band} ' BP Change %']); box off; axis square
        legend({'max./ amplifying phase','max./ suppressing phase'})
        ylim([-40 105])
    end
    
    subplot(2,3,3+((SScomb-1)*3))
    imagesc(AmpVec,EpsVec,squeeze(suppSave(1,:,:)))
    A= squeeze(suppSave(1,:,:));
    [~,ind] = max(A(:));
    [i,j] = ind2sub(size(A),ind); hold on
    scatter(AmpVec(j),EpsVec(i),100,'filled','Marker','o','MarkerEdgeColor','none','MarkerFaceColor','k')
    [~,ind] = min(A(:));
    [i,j] = ind2sub(size(A),ind); hold on
    scatter(AmpVec(j),EpsVec(i),100,'filled','Marker','v','MarkerEdgeColor','none','MarkerFaceColor','k')
    
    set(gca,'YDir','normal')
    xlabel('stim amplitude (relative)');
    ylabel('gating threshold (%)')
    title('B1 suppression')
    axis square; grid on; box off;
    colormap(brewermap(256,'*RdBu'))
    C(SScomb) = colorbar;
    if SScomb == 1
%         C(SSComb).Position = [0.920 0.605 0.009 0.314];
    else
%         C(SSComb).Position = [0.920 0.605 0.009 0.314]; 
    end
    caxis([-35 35])
end
    set(gcf,'Position',[ 290         235        1490         743])
    a = 1;

% ! shutdown /s

function bp = bpSpec(Hz,Pxy,flim)
bandInds = find(Hz>=flim(1) & Hz<= flim(2));
bp = sum(Pxy(bandInds))*(numel(bandInds)*(Hz(2)-Hz(1)));



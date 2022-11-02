function R = BAA_plotSweepSpectraBasic(R)
close all
rootan = fullfile(R.rootn, 'data', 'ConnectionSweep');

%% Set up plot options
R.CONnames = {'M2 -> STN','GPe -| STN'};
R.condname = {'Fitted','1% M2->STN','150% M2->STN','Fitted','1% STN->GPe','150% STN->GPe'};
cmap1 = brewermap(40,'Greys');
scmap = brewermap(4,'Set1');
statecmap{1} = scmap(1:2,:);
statecmap{2} = scmap(3:4,:);

% Loop through connections and plot spectra
figure(1)

for CON = 1:2
    load(fullfile(rootan, ['BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat.mat']),'feat')
    load(fullfile(rootan, ['BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1.mat']),'ck_1')
    
    % Plot M2 Spectra
    if CON == 1
        subplot(2,3,1)
    elseif CON == 2
        subplot(2,3,4)
    end
    plotSweepSpectra(R.frqz,feat,feat{1},cmap1,2:4:31,[1,1,1],statecmap{CON});
    title('M2 power')
    ylim([5e-5 1e-3]);%%         ylim([5e-17 5e-15])
    set(gca, 'YScale', 'log');
    
    % Plot STN Spectra
    if CON == 1
        subplot(2,3,2)
    elseif CON == 2
        subplot(2,3,5)
    end
    plotSweepSpectra(R.frqz,feat,feat{1},cmap1,2:4:31,[4,4,1],statecmap{CON})
    title('STN power')
    ylim([2e-4 1]);%     ylim([5e-16 1e-12])
    set(gca, 'YScale', 'log');
    
    % Plot STN/M2 Coherence
    if CON == 1
        subplot(2,3,3)
    elseif CON == 2
        subplot(2,3,6)
    end
    plotSweepSpectra(R.frqz,feat,feat{1},cmap1,2:4:31,[4,1,4],statecmap{CON})
    title('M2/STN coherence')
    ylim([0 1])
end
set(gcf,'Position',[ 518         250        1211         633])

%% Get spectral stats
for CON = 1:2
    load(fullfile(rootan, ['BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat.mat']),'feat')
    load(fullfile(rootan, ['BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1.mat']),'ck_1')
    load(fullfile(rootan, ['BB_' R.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim.mat']),'xsim')
    
    [~,~,bpowr,fpow,bcohr,fcoh,fpowCtx,bpowrCtx] = computeBetaSpectralStats(R.frqz,feat);
    ck_1 = ck_1(CON,:); % The scale for this connection modification
    ck_1 = ck_1(2:end);
    
    % Scale bpow to 0 (for plotting)
    [a zind] = min(abs(ck_1-1)); % base model
    bpowr = 100*(bpowr(2:end)-bpowr(1))/bpowr(1);
    bpowrCtx = 100*(bpowrCtx(2:end)-bpowrCtx(1))/bpowrCtx(1);
    bcohr = 100*(bcohr(2:end)-bcohr(1))/bcohr(1);
    fpow = fpow(2:end);
    fpowCtx = fpowCtx(2:end);
    fcoh = fcoh(2:end);
    
    powInds = find(bpowr>500);
    fpowCtx(powInds) = nan(1,numel(powInds));
    bpowrCtx(powInds) = nan(1,numel(powInds));
    fcoh(powInds) = nan(1,numel(powInds));
    bcohr(powInds) = nan(1,numel(powInds));
        
    %% Data Selection
    indsel = [1 30];
    dataSelect{CON} = xsim([1 2 31]);
    dataProperties(:,:,CON) = [bpowrCtx(indsel); fpowCtx(indsel); bpowr(indsel); fpow(indsel);  bcohr(indsel); fcoh(indsel);]
    
end
save(fullfile(rootan, ['BB_' R.out.tag '_DiscreteData.mat']),'dataSelect','dataProperties')


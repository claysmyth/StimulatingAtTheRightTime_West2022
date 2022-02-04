function [statBurstOverL sampBurstOverL] = burstOverlap(XF,fsamp,frqz,minper,band)

%% Define bursts at each population
XEnv = []; XPhi = [];
for ch = 1:6
    XEnv(ch,:) = abs(hilbert(XF(ch,:))); % retrieve amplitude envelope from filtered data
    % Set Threshold
    epsAmp(ch) = prctile(XEnv(ch,:),75); % get local eps for each population
    
    % Define Bursts
    ThreshX = double(XEnv(ch,:) > epsAmp(ch));
    minS = (minper/frqz(1))*fsamp; % set minimum to 3 periods of slowest
    betaBurstInds = SplitVec(find(ThreshX),'consecutive'); % Split up data based upon the target threshold
    segL = cellfun('length',betaBurstInds); % Find burst lengths
    burstSelInds = segL>minS; % Select bursts with above min length
    burstSelIndsSave{ch} = betaBurstInds(burstSelInds);
end

%% Now look for burst overlaps
burstOverl = nan(max(cellfun(@numel,burstSelIndsSave)),6,6);
sampBurstOverL = {};
for ch = 1:6
    for seg = 1:numel(burstSelIndsSave{ch})
        curBurst = burstSelIndsSave{ch}{seg};
        for ovch = setdiff(1:6,ch)
            ovBurst = [burstSelIndsSave{ovch}{:}];
            intInds = intersect(curBurst,ovBurst);
%             if ~isempty(intInds);
                burstOverl(seg,ovch,ch) = numel(intInds)./numel(curBurst);
                sampBurstOverL{ovch,ch}(seg) =  numel(intInds)./numel(curBurst);
%             else
%                 sampBurstOverL{ovch,ch}(seg) = nan;
%             end
        end
        
    end
end

% figure
% imagesc(squeeze(nanmean(burstOverl,1)),'AlphaData',~isnan(squeeze(nanmean(burstOverl,1))))
% caxis([0 1]); set(gca, 'ydir', 'reverse'); axis square;
% colorbar off
% 

statBurstOverL = nan(6,6,3);
if band == 1
    statBurstOverL(:,:,1) = triu(squeeze(nanmean(burstOverl,1)));
    statBurstOverL(:,:,2) = triu(squeeze(nanstd(burstOverl,[],1)));
    statBurstOverL(:,:,3) = triu(squeeze(sum(~isnan(burstOverl))));
    sampBurstOverL = (sampBurstOverL);
elseif band == 2
    statBurstOverL(:,:,1) = tril(squeeze(nanmean(burstOverl,1)));
    statBurstOverL(:,:,2) = tril(squeeze(nanstd(burstOverl,[],1)));
    statBurstOverL(:,:,3) = tril(squeeze(sum(~isnan(burstOverl))));
    
    sampBurstOverL = sampBurstOverL;
    
end


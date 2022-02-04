function initTime = relativePeakTiming(allEnv,twin,epsAmp,minper,frqz,fsamp)
minS = (minper/frqz(1))*fsamp; % set minimum to 3 periods of slowest
initTime = [];
for seg = 1:size(allEnv,3)
    Env = squeeze(allEnv(:,:,seg));
    for ch = 1:size(allEnv,1)
        burInd = SplitVec(find(Env(ch,:)>epsAmp(ch)),'consecutive'); % find all suprathreshold
            segL = cellfun('length',burInd); % criticise length for burst
        minlist = find(segL>minS,1,'first');
        if isempty(minlist)
            initTime(ch,seg) = nan;
        else
            initInd = burInd{minlist}(1)'; % time of first  burst
            initTime(ch,seg) = twin(initInd);
        end
    end
%     initTime(:,seg) = initTime(:,seg)-initTime(4,seg); % Normalize to STN onset
    
end

initTime ;
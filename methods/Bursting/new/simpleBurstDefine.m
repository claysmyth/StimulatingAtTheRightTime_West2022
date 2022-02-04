function [burstSelIndsout,XF,XEnv,XPhi,epsAmp,segLout,segAmpout] = simpleBurstDefine(X,fsamp,frqz,minper,eps,NSamp)
if nargin<5
    eps = 75;
end
% Filter Data
[dum,bpk,apk] = ft_preproc_bandpassfilter(X, fsamp,frqz,4,'but','twopass');
XF =filtfilt(bpk,apk,X')';


% Analytic Signal
XEnv = []; XPhi = [];
for ch = 1:6
    XEnv(ch,:) = abs(hilbert(XF(ch,:)));
    XPhi(ch,:) = angle(hilbert(XF(ch,:)));
    % Set Threshold
    epsAmp(ch) = prctile(XEnv(ch,:),eps);
end

for ch = 1:6
    % Define Bursts
    ThreshX = double(XEnv(ch,:) > epsAmp(ch));
    minS = (minper/frqz(1))*fsamp; % set minimum to 3 periods of slowest
    betaBurstInds = SplitVec(find(ThreshX),'consecutive'); % Split up data based upon the target threshold
    
    if eps == 0
        warning('Rewritting the burst definitions as eps = 0')
        for epoch = 1:NSamp(1)
            sind = randi(size(X,2)-NSamp(2),1);
            betaBurstInds{epoch} = sind:sind+NSamp(2)-1;
        end
    end
    
    segL = cellfun('length',betaBurstInds); % Find burst lengths
    burstSelInds = segL>minS; % Select bursts with above min length
    burstSelIndsout{ch} = betaBurstInds(burstSelInds);
    %
    
    segL = cellfun('length',burstSelIndsout{ch}); % Find burst lengths
    segLout{ch} = log10(segL);
    segAmp = cellfun(@(x) mean(XEnv(1:6,x),2),burstSelIndsout{ch},'UniformOutput',0);
    segAmpout{ch} = [segAmp{:}];
end

burstSelIndsout = burstSelIndsout{4};

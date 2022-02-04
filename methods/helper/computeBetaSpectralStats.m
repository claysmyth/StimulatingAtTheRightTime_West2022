function [bpowr_br,fpow_br,bpowr,fpow,bcohr,fcoh,fpowCtx,bpowrCtx] = computeBetaSpectralStats(Hz,feat)
bpow = []; fpow = [];
Hzlist = Hz(Hz>=14 & Hz<=30);
Hzlist_br = Hz(Hz>=14 & Hz<21);
for ck = 1:numel(feat)
    % Low Beta
    [dum b] = max(feat{ck}(1,4,4,3,Hz>=14 & Hz<21)); % Low Beta Power
    fpow_br(ck) = Hzlist_br(b);
    bpowr_br(ck) =  sum(feat{ck}(1,4,4,3,Hz>=14 & Hz<21))./sum(Hz>14 & Hz<21);
    
    % Whole Band
    [dum b] = max(feat{ck}(1,4,4,3,Hz>=14 & Hz<=30)); % Full Beta Power
    fpow(ck) = Hzlist(b);
    bpowr(ck) = sum(feat{ck}(1,4,4,3,Hz>=14 & Hz<=30))./sum(Hz>14 & Hz<30); % Full Beta Power
    
    
    % Coherence
    [peakcoh b] = max(feat{ck}(1,4,1,4,Hz>=14 & Hz<=30)); % Full Beta M2/STN Coh
    fcoh(ck) = Hzlist(b);
    bcohr(ck) = peakcoh; %A(b);
    
    % Whole Band M2
    [dum b] = max(feat{ck}(1,1,1,3,Hz>=14 & Hz<=30)); % Full Cortical Beta Power
    fpowCtx(ck) = Hzlist(b);
    bpowrCtx(ck) = sum(feat{ck}(1,1,1,3,Hz>=14 & Hz<=30))./sum(Hz>14 & Hz<30); % Full Beta Power
    
    
end
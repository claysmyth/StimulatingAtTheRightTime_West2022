function [twin,m2Env,stnEnv,dPhi,sw_twin,sw_PLV,outBurst,inBurst,allEnv,derv,intAmp]= getBurstTraces(burstSelInds,fsamp,XEnv,XPhi,Xpair)
oobInds = getRandOutofBurstInds(burstSelInds,size(XEnv,2)); % out of burst Inds
for seg = 1:numel(burstSelInds)
    zeroP = burstSelInds{seg}(1); % This is burst Onset
    Win = zeroP-(0.5*fsamp):zeroP+(0.75*fsamp);
    if Win(1)>1 && Win(end)<size(XEnv,2)
        derivWin = zeroP:zeroP+(0.5*fsamp);
        
        m2Env(:,seg) = XEnv(Xpair(1),Win);
        stnEnv(:,seg) = XEnv(Xpair(2),Win);
        allEnv(:,:,seg) = XEnv(:,Win);
        
        dPhi(:,seg) = wrapTo2Pi(diff(XPhi(Xpair,Win))); % M2/STN Phase
        sw_PLV(:,seg) = swPLV(dPhi(:,seg),0.15*fsamp); % Sliding Window PLV
        
        inBurst(seg) = circ_mean(wrapTo2Pi(diff(XPhi(Xpair,burstSelInds{seg})))');
        outBurst(:,seg) = circ_mean(wrapTo2Pi(diff(XPhi(Xpair,oobInds{seg})))');
        Yx = diff(unwrap(diff(XPhi(Xpair,derivWin))));
        derv(seg) = mean(abs(Yx)); % relative phase instability
        intAmp(seg) = sum(XEnv(Xpair(2),derivWin))/(numel(derivWin)*fsamp); %power in au/s-1
    end
end
twin = 1000.*linspace(-0.5,0.75,size(m2Env,1));
sw_twin = 1000.*linspace(-0.5,0.75,size(sw_PLV,1));

function oobInds = getRandOutofBurstInds(burstSelInds,N)
burstCatInds = [burstSelInds{:}];
for seg = 1:numel(burstSelInds)
    flag = 1;
    while flag
        sti = randi(N);
        candiInds = sti:sti+size(burstSelInds{seg},2)-1;
        
        % now check if out of burst
        if any(intersect(burstCatInds,candiInds)) || (candiInds(end)>N)
            flag = 1;
        else % accept
            oobInds{seg} = candiInds;
            flag = 0;
        end
    end
    
end



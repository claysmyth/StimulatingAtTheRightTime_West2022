function plotBurstTraces(twin,stnEnv,inBurst,dPhi,cmap)
% cmap = linspecer(3);
% cmap = brewermap(4,'RdYlBu');
oftab = [1.5 1.5 0.6 1; 2.5e-7 8e-7 0.6 1];
ostab = [1e-8 5e-8 0.05 0.01; 1e-8 2e-8 0.05 0.01];
lslist = {'-','--'};
npan = 4;
pvlist = nan(3,2,3);
for band = 1:2
    for state = 2:3   
        % Plot the mean within burst phase (radar plot
        subplot(2,npan,1+((band-1)*npan))
        magvec = 0.1 + (2.*0.025);
        % Radar plot
        Y = inBurst{band,state};
        X = inBurst{band,4};
        PS =circ_mean(Y,[],1);
        p(state) = plotPolarStem(PS,cmap(state,:)*1,magvec(1),'-','o',2,100);
        if state>1
            stattab(:,band,state) = radarStat(X,Y,cmap(state,:)*1,magvec(1),pi,'*',0.05);
        end
        title('Baseline Phase')
        rlim([0 0.16])
        
        % Plot the STN Burst Trace
        subplot(2,npan,2+((band-1)*npan))
        normvar = mean(stnEnv{band,state}(1,:)); %mean(stnEnv{band,1}(:));
        fy = mean(stnEnv{band,state}-normvar,2);
        fz = std(stnEnv{band,state}-normvar,[],2)./sqrt(size(stnEnv{band,state},2));
        boundedline(twin{band,state},fy,fz,lslist{band},'cmap',cmap(state,:),'alpha','transparency',0.5); hold on
        if state>1
            base = (stnEnv{band,4} - mean(stnEnv{band,4}(1,:)));
            fy = (stnEnv{band,state}-normvar);
            permSigStat(twin{band,state},base,fy,500,oftab(band,2),state*ostab(band,2),cmap(state,:),0.05,40)
        end
        title('STN Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); %ylim([-1 3]);
        grid on; box off; axis square
        
        % Plot the STN/M2 Relative Phase
        subplot(2,npan,3+((band-1)*npan))
        normvar = circ_mean(dPhi{band,state}(1,:),[],2);
        fy = circ_mean(dPhi{band,state},[],2)-normvar;
        fz = circ_std(dPhi{band,state},[],[],2)./sqrt(size(dPhi{band,state},2));
        boundedline(twin{band,state},wrapToPi(fy),fz,lslist{band},'cmap',cmap(state,:),'alpha','transparency',0.5); hold on
        if state>1
            base = wrapToPi(dPhi{band,4}-circ_mean(dPhi{band,4}(1,:),[],2));
            fy = wrapToPi(dPhi{band,state}-normvar);
            permSigStat(twin{band,state},base,fy,500,oftab(band,3),state*ostab(band,3),cmap(state,:),0.05,40)
        end
        title('M2/STN Phase'); ylabel('Relative Phase Angle'); xlabel('Time to Burst Onset (ms)'); %ylim([0.6 1.4]);
        grid on; box off; axis square
        
    end
end
l = legend(p(2:3),{'A','B'});
l.Position = [0.9005 0.0582 0.0463 0.0661];
        a = gca;
a.FontSize = 16;
a.FontName = 'Arial';

function stattab = radarStat(X,Y,cmap,magvec,shift,symb,alpha)
[pv tabww] = circ_wwtest(circ_mean(X,[],1),circ_mean(Y,[],1));
stattab = [circ_mean(X(:))-circ_mean(Y(:)) circ_std(X(:))-circ_std(Y(:)) tabww{4,2} tabww{2,5} tabww{2,6}];
if pv<alpha
    polarscatter(circ_mean(Y(:)',[],2)+shift,magvec,symb,'MarkerEdgeColor',cmap)
end

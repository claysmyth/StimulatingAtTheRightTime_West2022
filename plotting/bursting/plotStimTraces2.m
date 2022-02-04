function plotStimTraces2(twin,m2Env,stnEnv,sw_twin,sw_PLV,dPhi,middPhi)
cmap = brewermap(12,'Set1');
oftab = [4e-7 10e-7 1.2 1; 1e-8 6.5e-7 1.2 1];
ostab = [1e-8 5e-8 0.05 0.01; 1e-8 2e-8 0.05 0.01];
phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12

% cmap = brewermap(4,'RdYlBu');
philist = [100 6 12];
cmap = [0 0 0; cmap(philist(2)-1,:); cmap(philist(3)-1,:)];
lslist = {'-','--'};
npan = 5;

for band = 1:2
    for phi = 1:3
        if philist(phi) == 100
            stm = 1;
            pheff =  1;
        else
            stm = 2;
            pheff =  philist(phi);
        end
        subplot(2,npan,1+((band-1)*npan))
        normvar = mean(m2Env{band,pheff,stm}(1,:)); %mean(m2Env{band,1}(:));
        fy = mean(m2Env{band,pheff,stm}-normvar,2);
        fz = std(m2Env{band,pheff,stm}-normvar,[],2)./sqrt(size(m2Env{band,pheff,stm},2));
        boundedline(twin{band,pheff,stm},fy,fz,lslist{band},'cmap',cmap(phi,:),'alpha','transparency',0.5); hold on
        %         plot(twin{band,state},fy,'Color',cmap(state,:)); hold on;
        if phi>1
            base = (m2Env{band,1,1}-mean(m2Env{band,1,1}(1,:)));
            fy = (m2Env{band,pheff,stm}-normvar);
            permSigStat(twin{band,pheff,stm},base,fy,500,oftab(band,1),phi*ostab(band,1),cmap(phi,:),0.05,40)
        end
        title('M2 Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); % ylim([0 3]);
        grid on; box off; axis square
        
        subplot(2,npan,2+((band-1)*npan))
        normvar = mean(stnEnv{band,pheff,stm}(1,:)); %mean(stnEnv{band,1}(:));
        fy = mean(stnEnv{band,pheff,stm}-normvar,2);
        fz = std(stnEnv{band,pheff,stm}-normvar,[],2)./sqrt(size(stnEnv{band,pheff,stm},2));
        boundedline(twin{band,pheff,stm},fy,fz,lslist{band},'cmap',cmap(phi,:),'alpha','transparency',0.5); hold on
        if phi>1
            base = (stnEnv{band,1,1}-mean(stnEnv{band,1,1}(1,:)));
            fy = (stnEnv{band,pheff,stm}-normvar);
            permSigStat(twin{band,pheff,stm},base,fy,500,oftab(band,2),phi*ostab(band,2),cmap(phi,:),0.05,40)
        end
        title('STN Envelope'); ylabel('Normalized Amplitude'); xlabel('Time to Burst Onset (ms)'); %ylim([-1 3]);
        grid on; box off; axis square
        
        subplot(2,npan,3+((band-1)*npan))
        normvar = circ_mean(dPhi{band,pheff,stm}(1,:),[],2);
        fy = circ_mean(dPhi{band,pheff,stm},[],2)-normvar;
        fz = circ_std(dPhi{band,pheff,stm},[],[],2)./sqrt(size(dPhi{band,pheff,stm},2));
        boundedline(twin{band,pheff,stm},wrapToPi(fy),fz,lslist{band},'cmap',cmap(phi,:),'alpha','transparency',0.5); hold on
        if phi>1
            base = wrapToPi(dPhi{band,1,1}-circ_mean(dPhi{band,1}(1,:),[],2));
            fy = wrapToPi(dPhi{band,pheff,stm}-normvar);
            permSigStat(twin{band,pheff,stm},base,fy,500,oftab(band,3),phi*ostab(band,3),cmap(phi,:),0.05,0)
        end
        title('STN/M2 Phase'); ylabel('Relative Phase Angle'); xlabel('Time to Burst Onset (ms)'); %ylim([-1 3]);
        grid on; box off; axis square
        
        subplot(2,npan,4+((band-1)*npan))
        normvar = mean(abs(sw_PLV{band,pheff,stm}(1,:))); %mean(abs(sw_PLV{band,1}(:)));
        fy = mean(abs(sw_PLV{band,pheff,stm}),2);
        fz = std(abs(sw_PLV{band,pheff,stm}),[],2)./sqrt(size(sw_PLV{band,pheff,stm},2));
        boundedline(sw_twin{band,pheff,stm},fy,fz,'cmap',cmap(phi,:),'alpha','transparency',0.5); hold on
        if phi>1
            base = abs(sw_PLV{band,1,1});
            fy = abs(sw_PLV{band,pheff,stm});
            permSigStat(sw_twin{band,pheff,stm},base,fy,500,oftab(band,4),phi*ostab(band,4),cmap(phi,:),0.05,0)
        end
        
        title('M2/STN PLV'); ylabel('Normalized PLV'); xlabel('Time to Burst Onset (ms)'); %ylim([0.6 1.4]);
        grid on; box off; axis square
                
        % Radar plot
        subplot(2,npan,5+((band-1)*npan))
        magvec = 0.1 + ((phi-1).*0.025);
        
        % Radar plot
        Y = middPhi{band,pheff,stm};
        X = middPhi{band,1,1};
        PS =circ_mean(Y,[],1);
        p = plotPolarStem(PS,cmap(phi,:)*1,magvec(1),'-','o',3,100);
        if phi>1
%             stattab(:,band,phi) = radarStat(X,Y,cmap(phi,:)*1,magvec(1),pi,'*',0.05);
        end

        title('Baseline Phase')
        rlim([0 0.16])
        
    end
end
l = legend(p,{'Fitted','LB UP','LB Down'});
l.Position = [0.9005 0.0582 0.0463 0.0661];
function stattab = radarStat(X,Y,cmap,magvec,shift,symb,alpha)
[pv tabww] = circ_wwtest(circ_mean(X,[],1),circ_mean(Y,[],1));
stattab = [circ_mean(X(:))-circ_mean(Y(:)) circ_std(X(:)-Y(:)) tabww{4,2} tabww{2,5} tabww{2,6}];
if pv<alpha
    polarscatter(circ_mean(Y(:)',[],2)+shift,magvec,symb,'MarkerEdgeColor',cmap)
end

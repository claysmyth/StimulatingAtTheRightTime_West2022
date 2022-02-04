function a=  plotSweepSpectra(Hz,feat,featemp,cmap,condsel,chsel,stcmap)

% Mark out frequency borders of beta band
plot([14 14],[1e-32 1],'k--')
hold on
plot([21 21],[1e-32 1],'k--')
plot([30 30],[1e-32 1],'k--')


j = 0;
for i = condsel
    j = j+1;
    if max(squeeze(feat{i}(1,4,4,1,:)))<1e2; % If STN is over reasonable level
        a(i) = plot(Hz,squeeze(feat{i}(1,chsel(1),chsel(2),chsel(3),:)),'color',cmap(i,:),'LineWidth',2);
        hold on
    end
end
a(i+1) = plot(Hz,squeeze(featemp(1,chsel(1),chsel(2),chsel(3),:)),'k','LineWidth',2);
a(i+2) = plot(Hz,squeeze(feat{2}(1,chsel(1),chsel(2),chsel(3),:)),'Color',stcmap(1,:),'LineWidth',2);
a(i+3) = plot(Hz,squeeze(feat{31}(1,chsel(1),chsel(2),chsel(3),:)),'Color',stcmap(2,:),'LineWidth',2);

% legend(a(legsel),legn)
xlim([4 38])
xlabel('frequency (Hz)')
ylabel('amplitude (a.u.)')
title('Simulated STN Spectra')
box off
grid on
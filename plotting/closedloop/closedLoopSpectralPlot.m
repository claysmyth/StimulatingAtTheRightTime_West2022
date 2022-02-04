function closedLoopSpectralPlot(R,phaseShift,NcS,intpow,powspec_save,Hz)
cmap = brewermap(12,'Set1');
baseCon = find(NcS==1,1);

phaseShift = rad2deg(phaseShift);
for C =1:3
    if C == 1
        subplot(2,3,1) % Spectra Plot
        titbit = 'M2 Power';
    elseif C == 2
        subplot(2,3,2) % Spectra Plot
        titbit = 'STN Power';
    elseif C == 3
        subplot(2,3,3) % Spectra Plot
        titbit = 'STN/M2 Coherence';
    end
    phsel = 2:2:12;
    ip = 0;
    a = [];
    for i = phsel
        ip = ip+1;
        a(ip) = plot(Hz,squeeze(powspec_save(:,C,2,i,baseCon)),'color',cmap(i-1,:),'LineWidth',2);
        hold on
    end
    plot(Hz,squeeze(powspec_save(:,C,1,i,baseCon)),'color',[0 0 0],'LineWidth',2,'LineStyle','--');
    
    xlim([8 34])
%     legend(a,sprintfc('%.1f rad.',phaseShift(phsel)))
    xlabel('Frequency (Hz)'); ylabel([titbit])
    title([titbit])
    grid on; axis square
    
    % ARCs
    subplot(2,3,C+3) % Steady State Stats
    % Beta 1
    X = squeeze(intpow(C,1,:,:,baseCon))';
    X = 100.*(X(:,2)-X(:,1))./X(:,1);
    p(1) = plot(phaseShift,X,'color','k','LineWidth',1.5);
    hold on
    s = scatter(phaseShift(phsel),X(phsel),50,cmap(phsel-1,:),'filled');
    % Beta 2
    X = squeeze(intpow(C,2,:,:,baseCon))';
    X = 100.*(X(:,2)-X(:,1))./X(:,1);
    p(2) = plot(phaseShift,X,'color','k','LineWidth',1.5,'LineStyle','--');
    hold on
    s(2) = scatter(phaseShift(phsel),X(phsel),75,cmap(phsel-1,:),'filled');
    s(2).Marker = 'square';
    
    grid on;axis square
    
    xlabel('Stimulation Phase (radians)'); ylabel('Percentage Change')
    %     legend(p,{'\beta_1 (14-21 Hz) Power','\beta_2 (21-30 Hz) Power'})
    title([titbit ' Response Curve'])
    a = gca;
    a.XTick = rad2deg(([0 pi/2 pi 3*pi/2 2*pi]));
    
    xlim([0 360])
    ylim([-55 250]);
end
set(gcf,'Position',[ 711   418   957   560])

function plotConnectivityMats(R,connectMat,diffConnectCI,CON,cmap)
for state = 2:3
    if CON == 1
        sip = [0 1 3];
    else
        sip = [0 2 4];
    end
    A = abs(connectMat{1,state,CON})./abs(connectMat{1,1,CON});
    A(isnan(A)) = 0;
    B = abs(connectMat{2,state,CON})./abs(connectMat{2,1,CON});
    B(isnan(B)) = 0;
    AB = A+B;
    AB(logical(eye(size(AB)))) = nan;
    
    diffstat = squeeze(diffConnectCI{1,state,CON} + diffConnectCI{2,state,CON});
    AB(AB>squeeze(diffstat(1,:,:)) & AB<squeeze(diffstat(2,:,:))) = NaN;
    
    subplot(2,2,sip(state))
    imagesc(AB,'AlphaData',~isnan(AB))
    colormap(gca,cmap)
    a = gca;
    set(a, 'ydir', 'reverse');
    c = colorbar;
    ylabel(c, 'Change in PLV')
    if CON == 1
        caxis([0.60 1.40])
    elseif CON== 2
        caxis([0.4 1.6])
    end
    a.XTick = 1:6;
    a.XTickLabel =     R.chsim_name;
    a.XTickLabelRotation = -45;
    
    a.YTick = 1:6;
    a.YTickLabel =     R.chsim_name;
    a.XTickLabelRotation = -45;
    title(R.statename{CON,state})
    
    hold on
    plot([7 0],[7 0],'k','LineWidth',2)
    axis square;
    grid on
end
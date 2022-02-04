function plotStimExampleTrace(pU,xsim_ip,phangle,R)
for i = 1:3
    if i == 3
        cm = [0.2 0.2 0.2];
        fx{1} = zeros(1,numel(R.IntP.tvec_obs));
        fx{2} = xsim_ip{1,1}{1}{1}(1,1:end-1)';
        fx{3} = xsim_ip{1,1}{1}{1}(4,1:end-1)';
    elseif i == 2
        cm = R.cmap(i,:);
        fx{1} = pU{7,2}./50;
        fx{2} = xsim_ip{2,1}{phangle(i)}{1}(1,1:end-1)';
        fx{3} = xsim_ip{2,1}{phangle(i)}{1}(4,1:end-1)';
    elseif i == 1
        cm = R.cmap(i,:);
        fx{1} = pU{1,2}./50;
        fx{2} = xsim_ip{2,1}{phangle(i)}{1}(1,1:end-1)';
        fx{3} = xsim_ip{2,1}{phangle(i)}{1}(4,1:end-1)';
    end
    
       plot(R.IntP.tvec_obs,fx{1},'Color',cm,'LineWidth',2);
       hold on
       plot(R.IntP.tvec_obs,fx{2}-5e-6,'Color',cm,'LineWidth',2)
       plot(R.IntP.tvec_obs,fx{3}-10e-6,'Color',cm,'LineWidth',2)
end  
       xlim([9 10]); xlabel('Time (s)')
       a = gca;
       a.YTickLabel = [];
       a.YTick = [];
       box off
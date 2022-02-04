function BAA_plotStimTracesBase(R)
close all
rootan = [R.rootn 'data\phaseLockedStim'];

state = 1;
SScomb = 1;

CON = 1; % you only need the baseline model so either CON 1 or 2 are equal

load([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysisSweepState_' num2str(state) '_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],...
    'connectMat','specMat','statBurstOverMat','twin','m2Env','stnEnv','dPhi','sw_twin','sw_PLV',...
    'XS','pU','outBurst','inBurst')
% Plotting on base state only
Hz = R.frqz;
for stm = 1:2
    for phi = 1:12
        if stm == 1
            phiEf = 1;
        else
            phiEf = phi;
        end
        %                 spec = spectraSave{phiEf,stim};
        [F,Hz] = pwelch(XS{phiEf,stm}',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
        [C,Hz] = mscohere(XS{phiEf,stm}(1,:)',XS{phiEf,stm}(2,:)',1/R.IntP.dt,[],1/R.IntP.dt,1/R.IntP.dt);
        spec = [F C];
        
        powspec_save(:,:,stm,phi,1) = spec;
        intpow(:,1,stm,phi,1) = sum(spec(Hz>14 & Hz<=21,:))./(numel(find(Hz>21 & Hz<=30))*(Hz(2)-Hz(1)));
        maxpow(:,1,stm,phi,1) = max(spec(Hz>14 & Hz<=21,:));
        intpow(:,2,stm,phi,1) = sum(spec(Hz>21 & Hz<=30,:))./(numel(find(Hz>21 & Hz<=30))*(Hz(2)-Hz(1)));
        maxpow(:,2,stm,phi,1) = max(spec(Hz>21 & Hz<=30,:));
        intpow(:,3,stm,phi,1) = sum(spec(Hz>14 & Hz<=30,:))./(numel(find(Hz>21 & Hz<=30))*(Hz(2)-Hz(1)));
        maxpow(:,3,stm,phi,1) = max(spec(Hz>14 & Hz<=30,:));
        
    end
end

phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
phaseShift = phaseShift(1:12); %12

figure; closedLoopSpectralPlot(R,phaseShift,1,intpow,powspec_save,Hz)
figure; plotStimTraces2(twin,m2Env,stnEnv,sw_twin,sw_PLV,dPhi,inBurst)
cmap = brewermap(12,'Set1');
philist = [100 6 12];
cmap = [0 0 0; cmap(philist(2)-1,:); cmap(philist(3)-1,:)];
set(gcf,'Position',[401         321        1368         657])

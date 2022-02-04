function [R] = BAA_plotModelFit(Rorg)
close all
cmap = linspecer(6);
%% Load Model Parameters
load([Rorg.rootn 'data\modelfit\SimModelData_M10.mat'],'R','m','p')
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
% Simulate the Base Model
R.obs.trans.norm = 1; % Normalize the output spectra
R.obs.trans.gauss = 1; % Smooth
R.obs.obsstates = [1:6]; % All 6 nodes are observed
R.chloc_name = R.chsim_name; % Ensure sim names match to output names
R.obs.gainmeth = {'unitvar'};

%% Call the simulator
u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);
[r2,pnew,feat_sim,xsims,xsims_gl,wflag] = computeSimData120319(R,m,u,p,0);

% Plot Time Series
figure(1)
subplot(2,1,1)
X = xsims_gl{1};
for i = 1:6
    plot(R.IntP.tvec_obs,X(i,:)-(i*5),'Color',cmap(i,:))
    hold on  
end
xlim([31 33])

load([Rorg.rootn 'data\Storage\L6_lesion_rat_020317.mat'],'FTdata')

subplot(2,1,2)
X = FTdata.ContData.trial{1}([1 11 2 20],:);
T = FTdata.ContData.time{1};
for i = 1:4
    plot(T,X(i,:)-(i*0.2))
    hold on  
end
xlim([31 33])

% Plot Feature Space
R = prepareRatData_NoGauss_Group_NPD(R,0,0);

figure(2)
plotABCSpectraOnly(R.data.feat_xscale,R.data.feat_emp,feat_sim)
figure(3)
npdplotter_110717({R.data.feat_emp},{feat_sim},R.data.feat_xscale,R,[],[])
for i = 1:6; for j = 1:6; if i~=j;subplot(6,6,sub2ind([6 6],j,i));ylim([0 0.65]); end; end; end



function plotDataComparison(R);
% This function plots empirical data vs simulated
close all
load([ R.rootn '\data\Storage\L6_lesion_rat_020317.mat'])
ch = [1 12 5 19];
cmap = linspecer(4);

X = FTdata.ContData.trial{1}(ch,:);
for i = 1:4
    X(i,:) = normaliseV(X(i,:));
end
X = bandpass(X',[6 60],250)';

subplot(2,1,1)
p = plot(FTdata.ContData.time{1},X - ([0:4:12]'));
for l = 1:4
    p(l).Color = cmap(l,:);
end

xlim([20 25])

% Simulated Data
try
    load([ R.rootn '\data\rat_InDirect_ModelComp\ConnectionSweep\BB_InDrt_ModCompRev2_ConnectionSweep_CON_1_xsim_F1.mat'],'xsim')
catch
    warning('You need to simulate some data first!')
    return
end
% load([ R.rootn '\data\Storage\L6_lesion_rat_020317.mat'])
t = linspace(0,size(xsim{1}{1},2)/2000,size(xsim{1}{1},2));
X = xsim{1}{1}(1:4,:);
for i = 1:4
    X(i,:) = normaliseV(X(i,:));
end
X = bandpass(X',[6 60],1000)';

subplot(2,1,2)
p = plot(t,X - ([0:4:12]'));
for l = 1:4
    p(l).Color = cmap(l,:);
end

xlim([20 25])

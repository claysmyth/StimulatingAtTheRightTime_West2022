function [R,BB] = BAA_computeBurstWrapper(R)
rootan = fullfile(R.rootn, 'data', 'ConnectionSweep');
close all;

%% Sets up figure
% color options
scmap = brewermap(4,'Set1');
statecmap{1} = [0 0 0; scmap(1:2,:)];
statecmap{2} = [0 0 0; scmap(3:4,:)];

R.statename = {'Fitted','HD-Down','HD-Up'; 'Fitted','PS-Down','PS-Up'};

% Indices of pair to look at
Xpair = [1 4]; % [m2 stn]

segL = []; segA = [];

for CON = 1:2
    BB = [];
    load(fullfile(rootan, ['BB_' R.out.tag '_DiscreteData.mat']),'dataSelect')
    for state = 1:4
        %% Find Burst Epochs
        % Get Data
        if state<4
            X = dataSelect{CON}{state}{1};
        else
            X  = dataSelect{CON}{1}{1}; % state 4 is for permutations
        end
        for band = 1:2
            % Define Bands
            if band == 1 % B1
                butf = [14 21];
            elseif band == 2 %B2
                butf = [21 30];
            end
            
            % Set parameters
            fsamp = 1/R.IntP.dt;
            frqz = butf;
            minper = 3; % Min 3 cycles           
            [burstSelInds,XF,XEnv,XPhi,epsAmp,segL{band,state,CON},segA{band,state,CON}] = simpleBurstDefine(X,fsamp,frqz,minper);
            
            if state == 4 % This state is reserved for permutation (random out of burst)
                burstSelInds = permuteInds(burstSelInds,X);
            end
            
            % compute PLV
            [connectMat{band,state,CON},diffConnectCI{band,state,CON}] = computePLVConMat(XPhi,band);
            [statBurstPLVl{band,state,CON},~,] = burstPLV(XPhi,burstSelInds,band);
            
            % Compute Connectivity/Spectral Matrix
            [Pw pHz] = pwelch(X',fsamp,[],fsamp,fsamp);
            Pw = (Pw-mean(Pw(:)))./std(Pw(:));
            specMat{band,state,CON} = Pw(pHz>4 & pHz<=48,:);
            
            % Setup window in which to look
            winsize(1) = 0.075*fsamp; % pre onset
            winsize(2) = 0.075*fsamp; % post onset
            
            burstSelInds = cropBurstSelection(burstSelInds,winsize,size(X,2));% Remove bursts at ends
            
            [twin{band,state},m2Env{band,state},stnEnv{band,state},...
                dPhi{band,state},sw_twin{band,state},sw_PLV{band,state},...
                outBurst{band,state},inBurst{band,state},~,derv{band,state},intAmp{band,state}] = getBurstTraces(burstSelInds,fsamp,XEnv,XPhi,Xpair);
        end
    end
    
    % Plot traces
    figure(CON)
    plotBurstTraces(twin,stnEnv,inBurst,dPhi,statecmap{CON})
    
    % Plot Phase Instability vs Burst Power
    subplot(2,4,4)
    title('Phase Stability')
    cmap = statecmap{CON};
    for band = 1
        for state = [2 3]
            X = derv{band,state}; % phase instability
            Y = log10(intAmp{band,state}); % burst power
            scatter(X,Y,120,'MarkerFaceColor',cmap(state,:),'MarkerEdgeColor','none','MarkerFaceAlpha',0.3); hold on
            [r,p] = corrcoef(X,Y);
            if p(2)<0.05; % if correlation significant then plot regression
                [xCalc yCalc b Rsq] = linregress(X',Y');
                plot(xCalc,yCalc,'LineWidth',2,'Color',cmap(state,:));
            end
        end
    end
    xlabel('Average relative phase instability'); ylabel('Burst power'); grid on; box off
    axis square
    set(gcf,'Position',[ 0.0914    0.0722    1.3120    0.7104].*1e3)
    
    figure(102)
    if CON == 1
        cmap = brewermap(128,'RdBu');
    elseif CON == 2
        cmap = brewermap(128,'*PRGn');
    end
    plotConnectivityMats(R,statBurstPLVl,diffConnectCI,CON,cmap)
    
end
figure(102)
set(gcf,'Position',[401         321        1368         657])
save(fullfile(rootan, ['BB_' R.out.tag '_stateConnectivityMatrix.mat']), ...
    'specMat','connectMat','statBurstPLVl','segA','segL');

function burstSelIndsPerm = permuteInds(burstSelInds,X)
% This function makes a random out of burst draw
% it is  brute force, trial and error sampling
L = 1:size(X,2);
for i = 1:numel(burstSelInds)
    flag = 1;
    while flag
        p = randi(numel(L));
        pinds = p:p+numel(burstSelInds{i});
        if (pinds(end)+2000)<size(X,2) &&  (pinds(1)-2000)>1
            flag = 0;
        end
    end
    burstSelIndsPerm{i} = pinds;
end

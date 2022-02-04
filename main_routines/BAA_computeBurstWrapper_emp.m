function [R,BB] = BAA_computeBurstWrapper_emp(R)
rootan = 'C:\DATA\Rat_020317_processed';
% rootan = 'D:\Data\Rat_020317_processed';
close all;

R.statename = {'Fitted','HD-Down','HD-Up'; 'Fitted','PS-Down','PS-Up'};
anim{1} = {'C1','C2','C3','C4','C5','C6','C8'};
sfx{1} = '_control_rat_020317';
anim{2} = {'L4','L6','L13','L18','L19','L20','L22','L23'};
sfx{2} = '_lesion_rat_020317';

strlab = {'M1','STR','GPe','STN'};
for state =1:2
    BB = []; statBurstPLVl = {}; specMatAnim = {}; bmax = [];
    for aniN = 1:size(anim{state},2)
        load([rootan '\' anim{state}{aniN} sfx{state} '.mat'],'FTdata')
        fsamp = FTdata.fsample;
        X = FTdata.ContData.trial{1};
        X = resample(X',1/R.IntP.dt,fsamp)'; % resample to same as simulation
        fsamp = 1/R.IntP.dt; % new sample rate
        X = zscore(X,[],2);
        % find channel with max beta power
        Y = nan(6,size(X,2));
        for regName = 1:4
            chinds = find(strncmp(FTdata.ContData.label,strlab{regName},3));
            if numel(chinds)>0
                bbp = []; bbp = [];
                for i = 1:numel(chinds)
                    bbp(i) = bandpower(X(chinds(i),:),fsamp,[14 21]);
                end
                [bmax(aniN,regName),bind] = max(bbp);
                Y(regName,:) = X(chinds(bind),:);
            else
                Y(regName,:) = nan(1,size(X,2));
            end
        end
        X = Y;
        clear Y
        % Start feature computation
        for band = 1:2
            if band == 1
                fhz = 18; butf = [14 21];
            elseif band == 2
                fhz = 24; butf = [21 30];
            end
            
            % Set parameters
            frqz = butf;
            minper = 3; % Min 3 cycles
            
            [burstSelInds,XF,XEnv,XPhi,epsAmp] = simpleBurstDefine(X,fsamp,frqz,minper);
            PLVMat{band,aniN} = computePLVConMat(XPhi,band);
            
            [statBurstPLVl{band,aniN},~,] = burstPLV(XPhi,burstSelInds,band);
            
            [~,pHz] = pwelch(X(1:4,:)',fsamp,[],fsamp,fsamp);
            % Compute Connectivity/Spectral Matrix
            R = prepareRatData_NoGauss_Animal_NPD(R,aniN,state);
            Pww = [];
            for i = 1:6
                tmp = squeeze(R.data.feat_emp(1,i,i,1,:));
                Pww(:,i) = interp1(R.data.feat_xscale,tmp,pHz);
            end
            
            Pww = (Pww-nanmean(Pww(:)))./nanstd(Pww(:));
            specMatAnim{band,aniN} = Pww(pHz>4 & pHz<=48,:);
        end
    end
    %% Now average across animals
    for band = 1:2
        for aniN = 1:size(anim{state},2)
            connectMatEmp{band,state,1,aniN} = statBurstPLVl{band,aniN};
            connectMatEmp{band,state,1,aniN}(connectMatEmp{band,state,1,aniN}==1) = 0;
            
            plvMatEmp{band,state,1,aniN} = PLVMat{band,aniN}; 
            
            specMatEmp{band,state,1,aniN} = specMatAnim{band,aniN};
            
            % remove GPi/Thal
            plvMatEmp{band,state,1,aniN}(:,5:6) = nan;
            specMatEmp{band,state,1,aniN}(:,5:6) = nan;
            connectMatEmp{band,state,1,aniN}(:,5:6) = nan; connectMatEmp{band,state,1,aniN}(5:6,:) = nan;
        end
        
    end
end
rootan = [R.rootn 'data\ConnectionSweep'];
save([rootan '\BB_' R.out.tag '_EmpiricalStates_ConnectivityMatrix.mat'],'specMatEmp','connectMatEmp','plvMatEmp')

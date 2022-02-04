function [R,BB] = computeBurstWrapper_V3_empiricalData_perAnim(R)
% rootan = [R.rootn 'data\ConnectionSweep'];
% rootan = 'C:\DATA\Rat_020317_processed';
rootan = 'D:\Data\Rat_020317_processed';
close all;

R.statename = {'Fitted','HD-Down','HD-Up'; 'Fitted','PS-Down','PS-Up'};
anim{1} = {'C1','C2','C3','C4','C5','C6','C8'};
sfx{1} = '_control_rat_020317';
anim{2} = {'L4','L6','L13','L18','L19','L20','L22','L23'};
sfx{2} = '_lesion_rat_020317';

strlab = {'M1','STR','GPe','STN'};
for state =1:2
    BB = []; statBurstOverAnim = {}; statBurstPLVl = {}; connectMatAnim = {}; specMatAnim = {}; bmax = [];
    for aniN = 1:size(anim{state},2)
        load([rootan '\' anim{state}{aniN} sfx{state} '.mat'],'FTdata')
        fsamp = FTdata.fsample;
        X = FTdata.ContData.trial{1};
        X = resample(X',1/R.IntP.dt,fsamp)'; % resample to same as simulation
        fsamp = 1/R.IntP.dt; % new sample rate
        X = zscore(X,[],2);
        % find channel with max beta power
        Y = zeros(6,size(X,2));
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
            %             fsamp = 1/R.IntP.dt;
            frqz = butf;
            minper = 3; % Min 3 cycles
            %         peakHz = dataProperties(4,1,CON); % This loads in the peak STN frequency
            
            [burstSelInds,XF,XEnv,XPhi,epsAmp] = simpleBurstDefine(X,fsamp,frqz,minper);
            
            % compute burst overlaps
            [statBurstOverAnim{band,aniN}] = burstOverlap(XF,fsamp,frqz,minper,band);
            [statBurstPLVl{band,aniN},~,] = burstPLV(XPhi,burstSelInds,band);
            
            
            [connectMatAnim{band,aniN}] = computePLVConMat(XPhi,band);
            [~,pHz] = pwelch(X',fsamp,[],fsamp,fsamp);
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
    %      [~,aniSel] = max(bmax(:,4))
    for band = 1:2
        %         tmp1 = []; tmp2 = []; tmp3 = [];
        for aniN = 1:size(anim{state},2)
            %             tmp1(:,:,:,1) =  statBurstPLVl{band,aniN};% connectMatAnim{band,aniN};
            %             tmp2(:,:,1) =  specMatAnim{band,aniN};
            %             tmp3(:,:,:,1) =  statBurstOverAnim{band,aniN};
            connectMatEmp{band,state,1,aniN} = statBurstPLVl{band,aniN}; %mean(tmp1,4);
            connectMatEmp{band,state,1,aniN}(connectMatEmp{band,state,1,aniN}==1) = 0;
            specMatEmp{band,state,1,aniN} = specMatAnim{band,aniN};
            statBurstOverlEmp{band,state,1,aniN} = statBurstOverAnim{band,aniN};
            
            % remove GPi/Thal
            specMatEmp{band,state,1,aniN}(:,5:6) = nan;
            connectMatEmp{band,state,1,aniN}(:,5:6) = nan; connectMatEmp{band,state,1,aniN}(5:6,:) = nan;
        end
        
    end
end
rootan = [R.rootn 'data\ConnectionSweep'];
save([rootan '\BB_' R.out.tag '_EmpiricalStates_ConnectivityMatrix.mat'],'specMatEmp','connectMatEmp','statBurstOverlEmp')

function Pxy = logLogDetrend(F_scale,Pxy)
pwlinInds = (F_scale>48 & F_scale<52 );
Pxy(pwlinInds) = nan;
baseInds= (F_scale<6);
Pxy(baseInds) = nan;
tailinds = ((F_scale>6 & F_scale<=48 ) | (F_scale>52 & F_scale<98));
Pxy = log10(Pxy); F_scale = log10(F_scale);
[dum1 dum2 b Rsq] = linregress(F_scale(tailinds),Pxy(tailinds));
yCalc = [ones(length(F_scale),1) F_scale]*b;
Pxy = Pxy-yCalc;
Pxy = 10.^Pxy; F_scale = 10.^(F_scale);

truncinds = ((F_scale>98));
Pxy(truncinds) = nan;


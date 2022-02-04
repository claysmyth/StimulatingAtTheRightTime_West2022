function [R,BB] = BAA_computeStimAnalysis_sweep(R)
rootan = [R.rootn 'data\phaseLockedStim'];
% what type of stimulation
SScomb =1; % M2 stim phase locked to STN beta
lightweightFlag = 0; % reduces file size by not saving

Xpair = [1 4]; % The pair of pop indices to be analysed [m2 stn]

for CON = 1:2
            % Load in Data
        BB = [];
        load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_xsim' num2str(SScomb) '.mat'],'xsim_ip');
        %     load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_Rout' num2str(SScomb) '.mat'],'Rout');
        load([rootan '\BB_' R.out.tag '_phaseLockedStim_CON_' num2str(CON) '_pU_save' num2str(SScomb) '.mat'],'pU_save');
    
    for state = 1:31
        if state>1
            lightweightFlag = 1; % you dont need to save trajectories for state comparison
        end
        
        fsamp = 1/R.IntP.dt;
        for stm = 1:2
            % Set the phase range
            if stm==1
                philist = 1;
            elseif stm ==2
                philist = 1:12;
            end
            % Do base model
            Xbase = xsim_ip{1,state}{1};
            for phi = philist
                % Get Data
                X = xsim_ip{stm,state}{phi}{1};
                pU{phi,stm} = pU_save{state}{phi};
                
                % Get Burst Inds
                stimSelInds = SplitVec(find(abs(pU{phi,stm})>0),'consecutive'); % Split up data based upon stimulation gating
                
                % Setup window in which to look
                winsize(1) = 0.075*fsamp;
                winsize(2) = 0.075*fsamp; %0.15*fsamp;
                
                stimSelInds = cropBurstSelection(stimSelInds,winsize,size(X,2));% Remove bursts at ends
                
                XS{phi,1} = Xbase{1}(Xpair,[stimSelInds{:}]); % No stim
                XS{phi,2} = X(Xpair,[stimSelInds{:}]); % with stim
                
                
                for band = 1:2
                    if band == 1
                        fhz = 18; butf = [14 21];
                    elseif band == 2
                        fhz = 24; butf = [21 30];
                    end
                    
                    % Set parameters
                    fsamp = 1/R.IntP.dt;
                    frqz = butf;
                    minper = 3; % Min 3 cycles
                    
                    % Get Filtered/Analytic Signal
                    [~,XF,XEnv,XPhi,epsAmp] = simpleBurstDefine(X,fsamp,frqz,minper);
                    
                    % Get burst overlap
                    [statBurstOverMat{band,phi,stm}] = burstOverlap(XF,fsamp,frqz,minper,band);
                    
                    % Compute Connectivity Matrix
                    [Pw pHz] = pwelch(X',fsamp,[],fsamp,fsamp);
                    Pw = (Pw-mean(Pw(:)))./std(Pw(:));
                    specMat{band,phi,stm} = Pw(pHz>4 & pHz<=48,:);
                    connectMat{band,phi,stm} = burstPLV(XPhi,stimSelInds,band);
                    
                    [twin{band,phi,stm},m2Env{band,phi,stm},stnEnv{band,phi,stm},...
                        dPhi{band,phi,stm},sw_twin{band,phi,stm},sw_PLV{band,phi,stm},...
                        outBurst{band,phi,stm},inBurst{band,phi,stm},derv{band,phi,stm},intAmp{band,phi,stm}] = getBurstTraces(stimSelInds,fsamp,XEnv,XPhi,Xpair);
                    
                    disp([state band phi stm])
                end
            end
        end
        if ~lightweightFlag
            save([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysisSweepState_' num2str(state) '_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],...
                'connectMat','specMat','statBurstOverMat','twin','m2Env','stnEnv','dPhi','sw_twin','sw_PLV',...
                'XS','pU','outBurst','inBurst')
        else
            save([rootan '\BB_' R.out.tag '_phaseLockedStim_burstAnalysisSweepState_' num2str(state) '_CON_' num2str(CON) '_feat' num2str(SScomb) '.mat'],...
                'connectMat','specMat','statBurstOverMat');
        end
        
        disp(state)
    end
            clear xsim_ip

end



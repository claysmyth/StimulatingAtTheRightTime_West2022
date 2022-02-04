function [R] = BAA_sim_phaseLockedStim_SICompEps(Rorg)
close all
ip = 0;
load([Rorg.rootn 'data\modelfit\SimModelData_M10.mat'],'R','m','p')
R.rootn = Rorg.rootn;
R.filepathn = Rorg.filepathn;
warning('Loading Preloaded model, cant change simtime or model choice!!!')
for SScomb = 1:2
    EpsVec = [0:10:90];
    AmpVec = [1/8:1/8:1];
    for AmpN = 1:numel(AmpVec)
        for EpsN = 1:numel(EpsVec)
            if SScomb == 2
                %% Define stimulation conditions
                % Stimualating M2
                senssite = 1; % STN
                stimsite = 4; % M2
                stimFx = @zeroCrossingPhaseStim_v3;
                stim_sens = 'stimM2_sensSTN';
            elseif  SScomb == 1
                % Stimualating M2
                senssite = 4; % STN
                stimsite = 1; % M2
                stimFx = @zeroCrossingPhaseStim_v3;
                stim_sens = 'stimM2_sensSTN';
                
            end
            phflag = 1;
            % Setup stim parameterss
            R = typeIstimPars_v3(R);
            R.IntP.phaseStim.sensStm = [senssite stimsite];
            R.IntP.phaseStim.epsthresh = EpsVec(EpsN); % Overwrite
            R.IntP.phaseStim.stimGap = 0.2;
            R.IntP.phaseStim.stimperiod =1;
            R.IntP.phaseStim.stimAmp = AmpVec(AmpN);
            
            %% Connection Sets
            ck_1(1,:) = [1 logspace(-1,log10(5),34)];
            ck_1(2,:) = [1 logspace(-1,log10(1.90),34)];
                       
            % Simulation Coniditions
            R.obs.csd.df = 0.5;
            R = setSimTime(R,128);
            
            % Trans Options
            R.obs.trans.norm = 0;
            R.obs.gainmeth = {};
                        
            % Give all timeseries the same input - makes comparable
            rng(5453)
            uc = innovate_timeseries(R,m);
            XBase = p;
            
            % Phase To be Tested
            R.IntP.intFx = @spm_fx_compile_120319_stim;
            
            if phflag
                phaseShift = linspace(0,2.*pi,13); %13% List of phases to be tested
                phaseShift = phaseShift(1:12); %12
            else
                phaseShift = 0;
            end
            
            %% Setup Base Model
            Pbase = XBase;
            
            % Simulate Base Model
            uc_ip{1} = uc;
            R.IntP.phaseStim.switch = 0 ;
            R.IntP.phaseStim.phaseshift = 0;
            R.frqz = 2:0.2:150;
            R.IntP.compFx = @nullComp;
            [~,~,feat_sim_base{1},xsim_gl,xsim_ip_base{1}] = computeSimData120319(R,m,uc_ip{1},Pbase,0);
            % Work out the threshold
            R.IntP.phaseStim.eps = 0;
            [~,R] = zeroCrossingPhaseStim_v3([],R,0,xsim_gl{1},R.IntP.dt);
            
            %% Do Stimulation
            R.IntP.phaseStim.switch = 1;
            m = m; % initialise for parfor
            xsim_ip_stim = cell(1,12); feat_sim_stim = cell(1,12); pU = cell(1,12);
            %             parfor p = 1:numel(phaseShift)
            parfor p = 1:numel(phaseShift)
                Rpar = R;
                % Modulate the phase
                Rpar.IntP.phaseStim.phaseshift = phaseShift(p);
                Rpar.IntP.phaseStim.stimFx = stimFx;
                % Simulate with Stimulation
                [~,~,feat_sim_stim{p},~,xsim_ip_stim{p}]  = computeSimData120319(Rpar,m,uc_ip{1},Pbase,0);
                uexs = load([Rpar.rootn 'data\phaseStimSave\stim_tmp_' sprintf('%3.f',1000*Rpar.IntP.phaseStim.phaseshift)],'uexs');
                pU{p} = uexs.uexs(Rpar.IntP.phaseStim.sensStm(2),round(Rpar.obs.brn*(1/R.IntP.dt))+1:end);
                disp([AmpN EpsN p])
            end
            rmdir([R.rootn 'data\phaseStimSave\'],'s')
            % Create save version
            feat_sim_save{1,AmpN,EpsN} = feat_sim_base; feat_sim_save{2,AmpN,EpsN} = feat_sim_stim;
            xsim_ip{1,AmpN,EpsN} = xsim_ip_base; xsim_ip{2,AmpN,EpsN} = xsim_ip_stim;
            pU_save{AmpN,EpsN} = pU;
            
            ip = ip+1;
                       
            rootan = [Rorg.rootn 'data\phaseLockedStimEpsComp'];
            mkdir(rootan)
            
            save([rootan '\BB_' Rorg.out.tag '_phaseLockedStimEpsComp_EPS_feat' num2str(SScomb) '.mat'],'feat_sim_save','EpsVec','AmpVec')
        end
    end
end


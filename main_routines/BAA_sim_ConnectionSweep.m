function [R] = BAA_sim_ConnectionSweep(Rorg)

% Load in the saved model
load(fullfile(Rorg.rootn, 'data', 'modelfit', 'SimModelData_M10.mat'), ...
    'R','m','p')

% Sliding scale
% Used for: (1) Plot Spectra over sweeps
ck_1(1,:) = [1 logspace(log10(0.3),log10(3.11),30)]; % This is HD range
ck_1(2,:) = [1 logspace(log10(0.3),log10(1.33),30)]; % This is PS range

%% Trans Options
R.obs.SimOrd = 11;
R.obs.trans.norm = 0;
R.obs.gainmeth = {};

%% Compute Noise
% Give all timeseries the same input - makes comparable
uc = innovate_timeseries(R,m);
uc{1} = uc{1}.*sqrt(R.IntP.dt);
XBase = p;

%% Loop through Connections

% check if simulations exist

for CON = 1:2
    fprintf("\nCurrently running simulation for condition %d\n", CON);

    % create the output directory if does not exist
    rootan = fullfile(Rorg.rootn, 'data', 'ConnectionSweep');
    if ~exist(rootan, 'dir')
        mkdir(rootan)
    end

    % run simulation if saved outputs do not exist
    if ~exist(fullfile(rootan, ['BB_' Rorg.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat.mat']), "file")
        feat = {};
        xsim = {};
        parfor i = 1:size(ck_1,2)
            % Now modify parameter set
            Pbase = XBase;
            if CON == 1 % Hyperdirect
                Pbase.A{1}(4,1) = log(exp(Pbase.A{1}(4,1))*ck_1(CON,i)); %
            elseif CON == 2 % Pallidal-subthalamo
                Pbase.A{2}(4,3) = log(exp(Pbase.A{2}(4,3))*ck_1(CON,i)); %
            end
            [r2mean,pnew,feat_sim,dum1,xsim_gl] = computeSimData120319(R,m,uc,Pbase,0,1);
            feat{i} = feat_sim;
            xsim{i} = xsim_gl;

            disp([CON i])
        end


        save(fullfile(rootan, ['BB_' Rorg.out.tag '_ConnectionSweep_CON_' num2str(CON) '_feat.mat']),'feat')
        save(fullfile(rootan, ['BB_' Rorg.out.tag '_ConnectionSweep_CON_' num2str(CON) '_xsim.mat']),'xsim')
        save(fullfile(rootan, ['BB_' Rorg.out.tag '_ConnectionSweep_CON_' num2str(CON) '_ck_1.mat']),'ck_1')
    end
end


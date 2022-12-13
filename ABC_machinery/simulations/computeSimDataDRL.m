function [xsims, xsims_gl] = computeSimDataDRL(R_i, xsimsBuff_i, m_i, ...
    uc_i, uexs_i, p_i, tStepEnd_i, tvecRel_i)
% outputs
%   xsims - simulated dynamics (18, N)
%   xsims_gl - simulated local field potentials with right gain values

%% Simulate New Data
% Integrate in time master fx function
try
    [xsims, ~, wflag] = R_i.IntP.intFx(R_i, xsimsBuff_i, uc_i, uexs_i, ...
        p_i, m_i, tStepEnd_i);
catch
    disp('Simulation failed!')
    xsims{1} = nan(1,3);
    wflag = 1;
end

%% Run observation to get LFP

if wflag == 0
    try
        % Run Observer function
        % Subloop is local optimization of the observer gain
        glorg = p_i.obs.LF;
        gainlist = R_i.obs.glist;
        xsims_gl = cell(1,length(gainlist));
        for gl = 1:length(gainlist)
            p_i.obs.LF = glorg + gainlist(gl);
            if isfield(R_i.obs,'obsFx')
                [xsims_gl{gl}, ~, wflag(1)] = ...
                    R_i.obs.obsFx(xsims, m_i, p_i, R_i, tvecRel_i);
            else
                xsims_gl{gl} =xsims;
            end
            if any(wflag(1))
                error('Rejection at Observation!')
            end
        end
        xsims_gl = xsims_gl{1};

    catch
        disp('Observation/Cost Function Failure!')
        xsims_gl{1} = NaN;
    end
else
    disp('Sim Output contains NaNs!')
    xsims_gl{1} = NaN;
end

end
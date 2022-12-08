function [xsims_c R wflag] = DRL_observe_data(xstore,m,p,R, tvecRel)
    wflag = 0;

for condsel = 1:numel(R.condnames)
    xsims = xstore{condsel}(R.obs.outstates,:);
    % Delete burnin
    if size(xsims,2) > 1*round(R.obs.brn*(1/R.IntP.dt))
        if numel(tvecRel) > 0
            idxRemove = tvecRel < 0;
            xsims(:, idxRemove) = [];
        else
            warning('Debug mode: not removing brnIn period');
        end
    else
        wflag = 1;
        xsims_c{condsel} = xsims;
        warning('Simulation is shorter than burn length!!!')
        return
    end

    if any(isnan(xsims(:)))
        xsims_c{condsel} = xsims;
        warning('Simulation contains NaNs!!!')
        return
    end
    
    for i = 1:length(R.obs.gainmeth)
        switch R.obs.gainmeth{i}
            case 'obsnoise'
                CN = (R.obs.Cnoise.*exp(p.obs.Cnoise))';
                xsims = xsims + CN.*randn(size(xsims));
            case 'leadfield'
                LF = R.obs.LF.*exp(p.obs.LF);
                LFF = zeros(m.m);
                LFF(eye(size(LFF))~=0) = LF;
                xsims = LFF*xsims;
            case 'unitvar'
                for j = 1:size(xsims,1)
                    xsims(j,:) = (xsims(j,:) - R.DRL.chMean(j))./R.DRL.chStd(j);
                end
            case 'unitvarConcat'
                LM = [];
                for C = 1:numel(xstore)
                    LM = [LM xstore{C}(R.obs.outstates,round(R.obs.brn*(1/R.IntP.dt)):end)];
                end
                XM = mean(LM,2);
                XV = std(LM,[],2);
                xsims = (xsims-XM)./XV;            
        end
    end
    xsims_c{condsel} = xsims;
end

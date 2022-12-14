function X_comb = combineXStruct(vecX, R)

% loop through the codenames first
for condsel = 1:numel(R.condnames)
    % loop through the iterations
    for i = 1:numel(vecX)
        % if first loop then initialize
        if i == 1
            X_comb_cond = vecX{i}{condsel};
            continue
        end

        % otherwise pend to existing structure
        X_cond_curr = vecX{i}{condsel};
        
        % start with xsims
        fnXsims = fieldnames(X_comb_cond.xsims);
        for j = 1:numel(fnXsims)
            X_comb_cond.xsims.(fnXsims{j}) = [X_comb_cond.xsims.(fnXsims{j}), ...
                X_cond_curr.xsims.(fnXsims{j})];
        end

        % proceed to xsims_gl
        fnXsims_gl = fieldnames(X_comb_cond.xsims_gl);
        for j = 1:numel(fnXsims_gl)
            X_comb_cond.xsims_gl.(fnXsims_gl{j}) = [X_comb_cond.xsims_gl.(fnXsims_gl{j}), ...
                X_cond_curr.xsims_gl.(fnXsims_gl{j})];
        end

        % proceed to u
        fnU = fieldnames(X_comb_cond.u);
        for j = 1:numel(fnU)
            X_comb_cond.u.(fnU{j}) = [X_comb_cond.u.(fnU{j}), ...
                X_cond_curr.u.(fnU{j})];
        end

        % modify metadata
        X_comb_cond.metadata.numEpoch = X_comb_cond.metadata.numEpoch + ...
            X_cond_curr.metadata.numEpoch;

    end

    % append to output struct
    X_comb{condsel} = X_comb_cond;

end

end
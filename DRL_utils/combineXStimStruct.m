function X_stim_comb = combineXStimStruct(vecX_stim, R)

% loop through the codenames first
for condsel = 1:numel(R.condnames)
    % loop through the iterations
    for i = 1:numel(vecX_stim)
        % if first loop then initialize
        if i == 1
            X_stim_comb_cond = vecX_stim{i}{condsel};
            continue
        end

        % otherwise pend to existing structure
        X_stim_cond_curr = vecX_stim{i}{condsel};

        % proceed to xsims_gl
        fnXsims_stim_gl = fieldnames(X_stim_cond_curr.xsims_gl);
        for j = 1:numel(fnXsims_stim_gl)
            X_stim_comb_cond.xsims_gl.(fnXsims_stim_gl{j}) = ...
                [X_stim_comb_cond.xsims_gl.(fnXsims_stim_gl{j}), ...
                X_stim_cond_curr.xsims_gl.(fnXsims_stim_gl{j})];
        end

    end

    % append to output struct
    X_stim_comb{condsel} = X_stim_comb_cond;

end
end
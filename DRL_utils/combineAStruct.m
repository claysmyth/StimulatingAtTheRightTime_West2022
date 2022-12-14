function [A_comb] = combineAStruct(vecA, R)

% loop through the codenames first
for condsel = 1:numel(R.condnames)
    % loop through the iterations
    for i = 1:numel(vecA)
        % if first loop then initialize
        if i == 1
            A_comb_cond = vecA{i}{condsel};
            continue
        end

        % otherwise pend to existing structure
        A_cond_curr = vecA{i}{condsel};

        % proceed to xsims_gl
        fnUexs = fieldnames(A_cond_curr.uexs);
        for j = 1:numel(fnUexs)
            A_comb_cond.uexs.(fnUexs{j}) = [A_comb_cond.uexs.(fnUexs{j}), ...
                A_cond_curr.uexs.(fnUexs{j})];
        end

    end

    % append to output struct
    A_comb{condsel} = A_comb_cond;

end
end
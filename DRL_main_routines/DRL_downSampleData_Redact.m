function [X] = DRL_downSampleData_Redact(X, Rorg)

% make copy of input struct
R = deepCopyStruct(Rorg, 1);

% now downsample the X struct
dsFactor = fix(fix(1/R.IntP.dt) / R.DRL.dsFs);
for condsel = 1:numel(R.condnames)
    % unpack the structure
    X_cond = X{condsel};

    xsims_gl = X_cond.xsims_gl;

    % downsample the original state
    xsims_gl.S = cellfun(@(x) resample(x', 1, dsFactor)', xsims_gl.S, ...
        'UniformOutput', false);

    % append to outer structure
    X_cond.xsims_gl = xsims_gl;

    X{condsel} = X_cond;
end
end
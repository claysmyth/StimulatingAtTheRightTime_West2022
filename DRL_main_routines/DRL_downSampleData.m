function [X, X_stim, A] = DRL_downSampleData(X, ...
    X_stim, A, Rorg)

% make copy of input struct
R = deepCopyStruct(Rorg, 1);

% also downsample all data
dsFactor = fix(fix(1/R.IntP.dt) / R.DRL.dsFs);
for condsel = 1:numel(R.condnames)
    % unpack the structure
    X_cond = X{condsel};
    X_stim_cond = X_stim{condsel};
    A_cond = A{condsel};

    xsims_gl = X_cond.xsims_gl;
    xsims_gl_stim = X_stim_cond.xsims_gl;
    uexs = A_cond.uexs;

    % downsample the original state
    xsims_gl.S = cellfun(@(x) resample(x', 1, dsFactor)', xsims_gl.S, ...
        'UniformOutput', false);
    xsims_gl_stim.S = cellfun(@(x) resample(x', 1, dsFactor)', xsims_gl_stim.S, ...
        'UniformOutput', false);
    uexs.S = cellfun(@(x) resample(x, 1, dsFactor), uexs.S, ...
        'UniformOutput', false);

    % append to outer structure
    X_cond.xsims_gl = xsims_gl;
    X_stim_cond.xsims_gl = xsims_gl_stim;
    A_cond.uexs = uexs;

    X{condsel} = X_cond;
    X_stim{condsel} = X_stim_cond;
    A{condsel} = A_cond;
end
end
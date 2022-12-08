function vecStructOut = deepCopyStruct(structIn, n)
% make n deep copies of the same structure

% obtain list of fieldnames in input struct
fnStructIn = fieldnames(structIn);

% create n copies
for i = 1:n
    % copy all fields
    for j = 1:numel(fnStructIn)
        structOut_i.(fnStructIn{j}) = structIn.(fnStructIn{j});
    end

    % append to outer list
    vecStructOut{i} = structOut_i;
end

% check just one copy return struct directly
if n == 1
    assert(numel(vecStructOut) == 1);
    vecStructOut = vecStructOut{1};
end
end
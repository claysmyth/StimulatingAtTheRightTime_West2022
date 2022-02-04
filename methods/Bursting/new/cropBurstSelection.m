function burstSelInds = cropBurstSelection(burstSelInds,winsize,nmax)
    % Make sure the ends dont go off
    origin = cellfun(@(x) (x(1)-winsize(1)),burstSelInds);
    burstSelection1 = find(origin>1);
    
    origin = cellfun(@(x) (x(end)+winsize(2)),burstSelInds);
    burstSelection2 = find(origin<nmax);
    
    % Merge Selection
    burstSelection = intersect(burstSelection1,burstSelection2);
    burstSelection([1 end]) = []; % Remove first due to filter artefact
    burstSelInds = burstSelInds(burstSelection);
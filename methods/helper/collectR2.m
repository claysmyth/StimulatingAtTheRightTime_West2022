function R2 = collectR2(A,B)
for nf = 1:size(A,2)
    L(nf) = nansum((A{nf}(:)-B{nf}(:)).^2)./nansum((A{nf}(:)-nanmean(A{nf}(:))).^2);
end

R2 = 1-((sum(L(:)))/numel(L));
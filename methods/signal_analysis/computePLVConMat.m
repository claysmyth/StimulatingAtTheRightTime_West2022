function [conmat,diffCI] = computePLVConMat(XPhi,band)
cmblist = combnk(1:size(XPhi,1),2);
conmat = zeros(6);
diffCI = zeros(2,6,6);
for i=1:size(cmblist,1)
    XL = diff(XPhi([cmblist(i,1),cmblist(i,2)],:)); % relative phase
    
    winsize = floor(size(XL,2)/2000);
    r = rem(size(XL,2),2000);
    XL = reshape(XL(1:end-r),winsize,2000);
    
    
    if nargout>1
        % Do some statistics
        bootPLV = [];
        for cond = 1:2
            for lj = 1:500
                chi = randi(size(XPhi,2)-winsize);
                chj = randi(size(XPhi,2)-winsize);
                
                
                XPhi_i = XPhi(cmblist(i,1),chi:chi+winsize);
                XPhi_j = XPhi(cmblist(i,2),chj+winsize:-1:chj);
                
                bootPLV(lj,cond) = abs(PLV(XPhi_i-XPhi_j));
            end
        end
        
        bootdiff = 1+diff(bootPLV,1,2);
    else
        bootdiff = nan(100,1);
    end
    
    
    if band == 1 % rotate for each band so you can fit in one matrix
        conmat(cmblist(i,1),cmblist(i,2)) = mean(splitapply(@PLV,XL',1:size(XL,1)));
        diffCI(:,cmblist(i,1),cmblist(i,2)) = [prctile(bootdiff,1) prctile(bootdiff,99)];
    else
        conmat(cmblist(i,2),cmblist(i,1))  = mean(splitapply(@PLV,XL',1:size(XL,1)));
        diffCI(:,cmblist(i,2),cmblist(i,1)) = [prctile(bootdiff,1) prctile(bootdiff,99)];
    end
end
function [coeff,cmcoeff] = PLV(phi)
cmcoeff = sum(exp(1i*phi))./numel(phi);
% coeff = abs(cmcoeff);
coeff = (cmcoeff);
    
end
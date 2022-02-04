function [out,phiwin] = swPLV(phi,winsize)
slidePhi = slideWindow(phi,winsize,fix(winsize*0.95));
out = splitapply(@PLV,slidePhi,1:size(slidePhi,2));
end
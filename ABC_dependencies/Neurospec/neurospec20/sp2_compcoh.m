function [f,cl] = sp2_compcoh(sc1,cl1,sc2,cl2);
% function [f,cl] = sp2_compcoh(sc1,cl1,sc2,cl2)
% Function to implement difference of coherence test.
%
% Copyright (C) 2008, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%
% Based on Tanh-1(|Rab|) - Tanh-1(|Rcd|).
% Test requires same: samp_rate & seg_size.
% Also requires same number of segments: seg_tot.
% Test & confidence limit based on Rosenberg et al. (1989), P7.
% 
% Input arguments
%  sc1   Spectral coefficient matrix for  first coherence, Rab
%  cl1   cl structure for first coherence
%  sc2   Spectral coefficient matrix for second coherence, Rcd
%  cl2   cl structure for second coherence
%
% Output arguments
%  f    2 column matrix - 1: frequency, 2: Difference of coherence, Rab - Rcd.
%  cl   Structure with analysis parameters, including confidence limits.
%
% function [f,cl] = sp2_compcoh(sc1,cl1,sc2,cl2)

% Check numbers of arguments
if (nargin<4)
  error(' Not enough input arguments');
end  
if (nargout<2)
  error(' Not enough output arguments');
end

% Check estimates have identical Fourier frequencies.
if (cl1.seg_size~=cl2.seg_size)
  error(' Unequal length periodograms')
end
if (cl1.samp_rate~=cl2.samp_rate)
  error(' Unequal sampling rates')
end

% Check estimates have identical numbers of segmetns.
if (cl1.seg_tot_var~=cl2.seg_tot_var)
  error('Error - Unequal numbers of segments')
end

% Check sc matrices have correct dimension
if (size(sc1)~=[cl1.seg_size/2+1,3])
  error([' Error in size of matrix sc1'])
end  
if (size(sc2)~=[cl2.seg_size/2+1,3])
  error([' Error in size of matrix sc2'])
end  

% Generate output matrix
seg_size_2=(2:cl1.seg_size/2+1)';   % 2:N/2+1 Indexing for sc matrix (DC value not used).
% Column 1 - frequencies in Hz.
f(:,1)=(1:cl1.seg_size/2)'*cl1.df;
% Column 2 - Difference of coherence test.
rab=sc1(seg_size_2,3)./sqrt(sc1(seg_size_2,1).*sc1(seg_size_2,2));
rcd=sc2(seg_size_2,3)./sqrt(sc2(seg_size_2,1).*sc2(seg_size_2,2));
f(:,2)=atanh(abs(rab))-atanh(abs(rcd));

cl.seg_size=cl1.seg_size;         % Set values in cl structure. 
cl.seg_tot=cl1.seg_tot;           % Number of segments.
cl.seg_tot_var=cl1.seg_tot_var;   % Effective number of segments.
cl.samp_tot=cl1.samp_tot;         % Total number of samples.
cl.samp_rate=cl1.samp_rate;       % Sampling rate.
cl.df=cl1.df;                     % Deltaf.
% Confidence limit for comparison of coherence test.
cl.cmpcoh_c95=1.96*sqrt(1/cl.seg_tot_var);
cl.what='';

% Display No of segments & resolution 
disp(['Segments: ',num2str(cl.seg_tot),', Segment length: ',num2str(cl.seg_size)]);

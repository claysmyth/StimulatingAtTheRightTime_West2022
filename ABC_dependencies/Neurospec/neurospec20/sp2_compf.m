function [f,cl] = sp2_compf(sc1,cl1,ind1,sc2,cl2,ind2);
%function [f,cl] = sp2_compf(sc1,cl1,ind1,sc2,cl2,ind2);
% Function to implement log ratio difference of spectra test.
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
% Based on log_10 ratio with, supports different numbers of segments
% in original spectral estimates. Requires same samp_rate & seg_size
%
% Test & confidence limit based on Diggle (1990), p117-120.
% 
% Input arguments
%  sc1   Spectral coefficient matrix containing first spectrum, f11
%  cl1   cl structure for first spectrum
%  ind1  column in sc1 containing first spectral estimate (1 or 2)
%  sc2   Spectral coefficient matrix containing second spectrum, f22
%  cl2   cl structure for second spectrum
%  ind2  column in sc2 containing second spectral estimate (1 or 2)
%
% Output arguments
%  f    2 column matrix - 1: frequency, 2: log ratio difference of spectra (f11/f22)
%  cl   Structure with analysis parameters, including confidence limit.
%
%function [f,cl] = sp2_compf(sc1,cl1,ind1,sc2,cl2,ind2);

% Check numbers of arguments
if (nargin<6)
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

% Check sc matrices have correct dimension
if (size(sc1)~=[cl1.seg_size/2+1,3])
  error([' Error in size of matrix sc1'])
end  
if (size(sc2)~=[cl2.seg_size/2+1,3])
  error([' Error in size of matrix sc2'])
end  

% Check indices ind1, ind2
if (ind1<1 | ind1>2)
  error (' Input ind1 not valid')
end
if (ind2<1 | ind2>2)
  error (' Input ind2 not valid')
end

% Generate output matrix
seg_size_2=(2:cl1.seg_size/2+1)';   % 2:N/2+1 Indexing for sc matrix (DC value not used).
f(:,1)=(1:cl1.seg_size/2)'*cl1.df;                         % Column 1 - frequencies in Hz.
f(:,2)=log10(sc1(seg_size_2,ind1)./sc2(seg_size_2,ind2));  % Column 2 - Log ratio.

% Set values in cl structure.
cl.seg_size=cl1.seg_size;         % Segment length
cl.seg_tot_var1=cl1.seg_tot_var;  % Effective number of segments in f11
cl.seg_tot_var2=cl2.seg_tot_var;  % Effective number of segments in f22
cl.v1=2*cl.seg_tot_var1;          % Degrees of Freedom in f11
cl.v2=2*cl.seg_tot_var2;          % Degrees of Freedom in f22
cl.samp_tot=cl1.samp_tot;         % Total number of samples
cl.samp_rate=cl1.samp_rate;       % Sampling rate
cl.df=cl1.df;                     % Deltaf
% Upper and lower confidence limit for log ratio difference of spectra test.
cl.fcmp_c95u=log10(qf_v1v2(cl.v1,cl.v2,0.025));
cl.fcmp_c95l=-log10(qf_v1v2(cl.v2,cl.v1,0.025));
cl.what='';

% Display DOF and segment length.
disp(['v1 v2: ',num2str(cl.v1),', ',num2str(cl.v2),', Segment length: ',num2str(cl.seg_size)]);

return

function f = qf_v1v2(v1,v2,p)

% function f = qf_v1v2(v1,v2,p)
%
% function to estimate value of F distribution which satisfies:
% Q(F|v1,v2)=p, for specified v1 and v2.
% Fixed values of p supported: 0.01, 0.025, 0.05 (1%, 2%, 5%).
% Page references are to Abramowitz & Stegun, Ch 26.

if (v1<2 | v2<2)
  error([' Error in function qf_v1v2 - Incorrect v1, v2: ',num2str(v1),', ',num2str(v2)]);
end  

% See Example 24, p961
a=v2/2;
b=v1/2;

% h, lambda and w part of 26.5.22, p945.
h=2/(1/(2*a-1)+1/(2*b-1));
switch (p)              % y for fixed percentage points of normal dist. 
  case 0.01             % 1.0 % level
    y=2.32635;
  case 0.025            % 2.5 % level
    y=1.95996;
  case 0.05             % 5.0 % level
    y=1.64485;
  otherwise
    error ([' Error in function qf_v1v2 - Value of p not supported: ',num2str(p)]);
end
lambda=(y*y-3)/6;
w=(y*sqrt(h+lambda)/h) - (1/(2*b-1)-1/(2*a-1))*(lambda+5/6-2/(3*h));

% f estimated from 26.6.16, p947.
f=exp(2*w);

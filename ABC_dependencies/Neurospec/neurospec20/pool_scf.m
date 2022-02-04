function [out_f,out_v] = pool_scf(sc,cl,in_f,in_v);
% [out_f,out_v] = pool_scf(sc,cl,in_f,in_v);
% Function to esimtate pooled spectral coefficients and calculate
%  chi-squared difference tests on spectra and coherence.
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
% Input arguments
%  sc    spectral coefficient matrix, from sp2* analysis
%  cl    cl structure from sp2* analysis
%  in_f  Input pooled spectral matrix (Optional - omit for new pooled analysis). 
%  in_v  Input pooled spectral variable structure (Optional - omit for new pooled analysis).
%
% Output arguments
%  out_f  Updated pooled spectral matrix (see comments in code for contents).
%  out_v  Updated pooled spectral variable structure (see comments in code for contents).
%
% Reference for pooled spectral analysis:
% Amjad, A.M., Halliday, D.M., Rosenberg J.R. & Conway B.A (1997).
%  An extended difference of coherence test for comparing and combining several
%   independent coherence estimates - theory and application to the study of
%   motor units and physiological tremor.
%  Journal of Neuroscience Methods, 73(1), 69-79.
%
% [out_f,out_v] = pool_scf(sc,cl,in_f,in_v);

% Note: JNM refers to above Journal of Neuroscience Methods article.

% Check numbers of arguments. Both in_f and in_v should be specified. 
if (nargin<2 | nargin==3)
  error(' Not enough input arguments');
end  
if (nargout<2)
  error(' Not enough output arguments');
end

% Check sc and cl have correct dimensions
if (length(size(sc))~=2)
  error(['Matrix sc - incorrect number of dimensions'])
end

if (length(cl)~=1)
  error(['Structure cl - incorrect number of dimensions'])
end

% Check sizes of input matrices are compatible with pooled analysis.
if (size(sc)~=[cl.seg_size/2+1,3])
  error(['Error in size of matrix sc'])
end  

% New pooled analysis: No pooled variables input  or  in_f is empty matrix.
% Set up matrix and structure for new pooled coherence analysis
if (nargin<3 | (nargin==4 & isempty(in_f)))
  in_f=zeros(cl.seg_size/2,13);  % Has seg_size/2 rows - zero frequency component is not pooled.
  in_v.weight_tot=0;             % Sum of weightings used, Sum(L_i) in JNM. 
  in_v.fil_tot=0;                % Number of files included.
  in_v.seg_tot=0;                % Number of segments of data pooled. 
  in_v.samp_tot=0;               % Number of samples of data pooled.
  in_v.samp_rate=cl.samp_rate;   % Sampling rate  (Should be same for each contributing file).
  in_v.seg_size=cl.seg_size;     % Segment length (Should be same for each contributing file).
end  

if (size(in_f)~=[cl.seg_size/2,13])
  error(['Error in size of matrix in_f'])
end  

% Check analysis parameters against previous files for identical Fourier frequencies.
if (cl.seg_size~=in_v.seg_size)
  error(['Error - Unequal length periodograms: ',num2str([cl.seg_size in_v.seg_size])])
end
if (cl.samp_rate~=in_v.samp_rate)
  error(['Error - Unequal sampling rates: ',num2str([cl.samp_rate in_v.samp_rate])])
end

% Indexing for sc matrix not including component at zero frequency.
seg_size_2=(2:cl.seg_size/2+1)';

% Estimate: sqrt(f11*f22), magnitude of coherency, coherence and z transform.
f_fac=sqrt(sc(seg_size_2,1).*sc(seg_size_2,2));
chy=abs(sc(seg_size_2,3))./f_fac;
ch=chy.*chy;
z=atanh(chy);

% Update histogram counting significant coherence values.
out_f(:,1)=in_f(:,1);
ch_ind=find(ch>cl.ch_c95);         % Include all values above 95% confidence limit.
out_f(ch_ind,1)=out_f(ch_ind,1)+1; % Column 1 is histogram count of significant coherences.

% Weighting factor is effective number of segments (L_i in JNM)
pool_fac=cl.seg_tot_var;

% Log_10 transformed spectra: log10(f11), log10(f22)
f11=log10(sc(seg_size_2,1));
f22=log10(sc(seg_size_2,2));

% Update pooled spectral coefficient matrix.                     % Col  1 is count of significant coherence.
out_f(:,2) =in_f(:,2) +sc(seg_size_2,1)*pool_fac;                % Col  2 is pooled f11       (2.11) in JNM.
out_f(:,3) =in_f(:,3) +sc(seg_size_2,2)*pool_fac;                % Col  3 is pooled f22       (2.11) in JNM.
out_f(:,4) =in_f(:,4) +real(sc(seg_size_2,3))*pool_fac;          % Col  4 is pooled Re{f21}   (2.11) in JNM.
out_f(:,5) =in_f(:,5) +imag(sc(seg_size_2,3))*pool_fac;          % Col  5 is pooled Im{f21}   (2.11) in JNM.
out_f(:,6) =in_f(:,6) +z*pool_fac;                               % Col  6 is Sum(L_i * z_i)   (2.10) in JNM.
out_f(:,7) =in_f(:,7) +z.*z*pool_fac;                            % Col  7 is Sum(L_i * z_i^2) (2.10) in JNM.
out_f(:,8) =in_f(:,8) +(pool_fac*real(sc(seg_size_2,3))./f_fac); % Col  8 is pooled Re{r21}
out_f(:,9) =in_f(:,9) +(pool_fac*imag(sc(seg_size_2,3))./f_fac); % Col  9 is pooled Im{r21}
out_f(:,10)=in_f(:,10)+f11*pool_fac;                             % Col 10 is Sum(L_i * f11_i)
out_f(:,11)=in_f(:,11)+f11.*f11*pool_fac;                        % Col 11 is Sum(L_i * f11_i^2)
out_f(:,12)=in_f(:,12)+f22*pool_fac;                             % Col 12 is Sum(L_i * f22_i)
out_f(:,13)=in_f(:,13)+f22.*f22*pool_fac;                        % Col 13 is Sum(L_i * f22_i^2)

% Update pooled spectral variable structure
out_v.weight_tot=in_v.weight_tot+pool_fac;      % Total number of effective segments
out_v.fil_tot=in_v.fil_tot+1;                   % Total number of files analysed
out_v.seg_tot=in_v.seg_tot+cl.seg_tot;          % Total number of segments analysed
out_v.samp_tot=in_v.samp_tot+cl.samp_tot;       % Total number of data samples analysed
out_v.samp_rate=cl.samp_rate;                   % Sampling rate
out_v.seg_size=cl.seg_size;                     % Segment length

% Display Running count: No of files & segments
disp(['Pooled - Files: ',num2str(out_v.fil_tot),', Segments: ',num2str(out_v.seg_tot)]);

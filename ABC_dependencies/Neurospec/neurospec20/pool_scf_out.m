function [f,t,cl,sc] = pool_scf_out(plf,plv);
% [f,t,cl,sc] = pool_scf_out(plf,plv);
% Function to convert pooled spectral coefficients into f, t, cl for plotting.
% Input arguments plf, plv are output of function pool_scf.
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
%   plf   Pooled spectral matrix (output from pool_scf).
%   plv   Variables from pooled analysis (output from pool_scf).
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%  sc   optional matrix for spectral coefficients.
%
% Output parameters
%  f column  1  frequency in Hz.
%  f column  2  Log pooled (input)  spectrum.
%  f column  3  Log pooled (output) spectrum.
%  f column  4  Pooled coherence (calculated from pooled coherency).
%  f column  5  Pooled phase.
%  f column  6  Pooled coherence (calculated from pooled spectra).
%  f column  7  Coherence histogram.
%  f column  8  Chi-squared difference of coherence test.
%  f column  9  Chi-squared difference of  input spectrum test.
%  f column 10  Chi-squared difference of output spectrum test.
%
%  t column 1  Lag in ms.
%  t column 2  Pooled cumulant density.
%
%  cl.seg_size    Segment length.
%  cl.seg_tot     Number of segments.
%  cl.seg_tot_var Effective number of segments, used to calculate confidence limits.
%  cl.samp_tot    Number of samples analysed.
%  cl.samp_rate   Sampling rate of data (samps/sec).
%  pl_fil_tot     Number of files.
%  pl_weight_tot  Sum of weightings used, Sum(L_i).
%  cl.dt          Time domain bin width (ms).
%  cl.df          Frequency domain bin width (Hz).
%  cl.f_c95       95% confidence limit for Log spectral estimates.
%  cl.ch_c95      95% confidence limit for coherence.
%  q_c95          95% confidence limits for cumulant.
%  chi_c95        95% confidence limit for Chi^2 test, JNM (2.9)
%
% [f,t,cl,sc] = pool_scf_out(plf,plv);

% Definition of upper 95% confidence limit for chi^2 distribution with DOF 1 - 30.
chi_005=[3.84,5.99,7.81,9.49,11.1,12.6,14.1,15.5,16.9,18.3,19.7,21,22.4,23.7,25, ...
         26.3,27.6,28.9,30.1,31.4,32.7,33.9,35.2,36.4,37.7,38.9,40.1,41.3,42.6,43.8];

% Variance factor for chi-squared comparison of spectra test (~5.302)
a=1/(log10(exp(1)))^2;

% Normalise pooled coefficients in columns 2:5, 8:9 by plv.weight_tot
plf(:,[2:5 8 9])=plf(:,[2:5 8 9])./plv.weight_tot;

% Pooled coherence based on pooled spectra
coh_f=((plf(:,4).*plf(:,4))+(plf(:,5).*plf(:,5)))./(plf(:,2).*plf(:,3));

% Pooled coherence based on pooled coherency
coh_r=(plf(:,8).*plf(:,8))+(plf(:,9).*plf(:,9));

deltaf=plv.samp_rate/plv.seg_size;        % Resolution - spacing of Fourier frequencies in Hz.

% Construct output spectral matrix f.
seg_size_2=(1:plv.seg_size/2)';           % Indexing for output, DC component not output.
f(:,1)=(seg_size_2)*deltaf;               % Column 1 - frequencies in Hz.
f(:,2)=log10(plf(:,2));                   % Column 2 - Log spectrum ch 1.
f(:,3)=log10(plf(:,3));                   % Column 3 - Log spectrum ch 2.
f(:,4)=coh_r;                             % Column 4 - Coherence (Pooled coherency estimate).
f(:,5)=angle(complex(plf(:,4),plf(:,5))); % Column 5 - Phase.
f(:,6)=coh_f;                             % Column 6 - Coherence (Pooled spectra estimate).
f(:,7)=plf(:,1)./plv.fil_tot;             % Column 7 - Coherence histogram (Normalised).
f(:,8) =2*(plf(:,7) -(plf(:,6) .*plf(:,6) ./plv.weight_tot)); % Column  8 chi^2 test coherence.
f(:,9) =a*(plf(:,11)-(plf(:,10).*plf(:,10)./plv.weight_tot)); % Column  9 chi^2 test f11.
f(:,10)=a*(plf(:,13)-(plf(:,12).*plf(:,12)./plv.weight_tot)); % Column 10 chi^2 test f22.

% Estimate cumulant density using inverse DFT of cross spectrum.
deltat=1000.0/plv.samp_rate; % dt in msec.
seg_size=plv.seg_size;       % Segment length in DFT: T.

% Construct complex valued matrix to use in ifft routine
f21(1,1)=0+0i;                                    % First value is zero (DC component)
f21(2:seg_size/2+1,1)=complex(plf(:,4),plf(:,5)); % Next T/2 values copied straight over.
% Index in reverse order (Frequency values: T/2-1 to 1; Indexing T/2 to 2)
tmp_ind=(seg_size/2:-1:2)';
% Copy remaining T/2-1 values in reverse order, conjugate Im part.
f21(seg_size/2+2:seg_size,1)=complex(plf(tmp_ind,4),-plf(tmp_ind,5));

cov=ifft(f21);              % Inverse DFT.

% Construct output time domain matrix t.
% Column 1 - time in msec. Range (-T/2)*dt to (T/2-1)*dt.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;

% Column 2 - cumulant, shifted by T/2 so that time zero is in centre.
% NB 2pi/T factor is 2*pi, since ifft includes 1/T term.
t([seg_size/2+1:seg_size,1:seg_size/2],2)=real(cov(1:seg_size))*2*pi; % PBMB (5.9).

% Estimate variance of cumulant density estimate.
var_fac=4*pi*pi/(seg_size*plv.samp_tot);                       % Factor (2pi/T)(2pi/R).
q_var=var_fac*2*sum(plf(1:seg_size/2,2).*plf(1:seg_size/2,3)); % PBMB (6.10), No DC value.

% Confidence limits for chi^2 difference of coherence test.
chi_df=plv.fil_tot-1;        % DOF in chi^2 variate.
if (chi_df<length(chi_005))
  chi_c95=chi_005(chi_df);   % Value in above look up table, for DOF <= 30.
else
  t1=2.0/(9.0*chi_df);       % For DOF > 30, estimate directly.
  t2=1.0-t1+1.6449*sqrt(t1); % Calculated using cube approximation,
  t3=t2*t2*t2;               % Ref: Abramowitz & Stegun (26.4.17), P985.
  chi_c95=chi_df*t3;
end

% Construct optional spectral matrix sc, DC value of 0 is pre appended to output for DC component.
if (nargout>3)
  sc(:,1)=[0;plf(:,2)];                        % Column 1 - f11
  sc(:,2)=[0;plf(:,3)];                        % Column 2 - f22
  sc(:,3)=complex([0;plf(:,4)],[0;plf(:,5)]);  % Column 3 - f21 (Complex valued)
end

% Construct cl structure, includes confidence limits for parameter estimates.
cl.seg_size=plv.seg_size;         % Set values in cl structure. 
cl.seg_tot=plv.seg_tot;           % No of segments.
cl.seg_tot_var=plv.weight_tot;    % Effective no of segments.
cl.samp_tot=plv.samp_tot;         % No of samples.
cl.samp_rate=plv.samp_rate;       % Sampling rate.
cl.pl_fil_tot=plv.fil_tot;        % No of files.
cl.pl_weight_tot=plv.weight_tot;  % Sum of weightings used, Sum(L_i).
cl.dt=deltat;                     % dt.
cl.df=deltaf;                     % df.
cl.f_c95=0.8512*sqrt(1/plv.weight_tot);    % Confidence limit for spectral estimates, PBMB (6.2).
% N.B. Confidence interval for log plot of spectra is TWICE this value.
cl.ch_c95=1-0.05^(1/(plv.weight_tot-1));   % Confidence limit for coherence, PBMB (6.6).
cl.q_c95=1.96*sqrt(q_var);                 % Confidence limits for cumulant, PBMB (6.11).
cl.chi_c95=chi_c95;                        % Confidence limit for Chi^2 test, JNM (2.9)

cl.N1=0;                          % N1, No of events in ch 1. (set to zero for pooled analysis)
cl.N2=0;                          % N2, No of events in ch 2.           "
cl.P1=0;                          % P1, mean intensity ch1.             "
cl.P2=0;                          % P2, mean intensity ch2.             "
cl.opt_str='Pooled';              % Set options string to indicate pooled analysis.
cl.what='';                       % Field for plot label.

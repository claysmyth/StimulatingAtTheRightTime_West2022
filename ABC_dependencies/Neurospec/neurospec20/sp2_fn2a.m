function [f,t,cl,sc] = sp2_fn2a(d1,d2,seg_samp_no,samp_rate,flags);
% function [f,t,cl,sc] = sp2_fn2a(d1,d2,seg_samp_no,samp_rate,flags);
% Function with core routines for periodogram based spectral estimates.
% Implements a weighted periodogram analysis with:
%  Variable number of data points in each segment, T
%  Fixed DFT segment length, S.
% see NOTE below about calling this routine.
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
%  Inputs are two matrices containing pre processed time series or point process data.
%  Matrices have L columns with S rows.
%   L = Number of segments: seg_tot.
%   S = DFT segment length: seg_size.
% 
% Input arguments
%  d1          Channel 1 data matrix.
%  d2          Channel 2 data matrix.
%  seg_samp_no Vector with number of data points (T) in each segment.
%  samp_rate   Sampling rate (samples/sec)
%  flags       Structure with flags to control processing options.
%              Flags supported:
%                line: Apply simple filter to suppress mains/line frequency (0:No; 1:Yes).
%                inv:  Invert channel reference for phase & cumulant        (0:No; 1:Yes).
%                han:  Additional smoothing to spectra using hanning filter (0:No; 1:Yes).
%                wind: Apply data window (split cosine taper) to raw data.
%                        (0:No; 1,2,3,...: Type).  1:10%; 2:20%; 3:50%; 4:100% taper.
%                gain: Estimate additional system identification parameters (0:No; 1:Yes).
%                        Log{|Gain|} included in f matrix.
%                        Impulse response included in t matrix.
%                cov:  Calculates Periodogram COV test for both channels.
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%  sc   optional matrix for spectral coefficients.
%
% Output parameters
%  f column 1       frequency in Hz.
%  f column 2       Log  input/d1 spectrum.
%  f column 3       Log output/d2 spectrum.
%  f column 4       Coherence.
%  f column 5       Phase.
%  f column 6       Log of gain magnitude (with s option).
%  f column 6 or 7  Periodogram COV test  input channel (with c option).
%  f column 7 or 8  Periodogram COV test output channel (with c option).
%
%  t column 1       Lag in ms.
%  t column 2       Cumulant density.
%  t column 3       Impulse response (with s option).
%
%  cl.type          Analysis type (0, 1, 2)
%  cl.seg_size      Segment length.
%  cl.seg_tot       Number of segments.
%  cl.seg_tot_var   Effective no of segments, used to calculate confidence limits.
%  cl.samp_tot      Number of samples analysed.
%  cl.samp_rate     Sampling rate of data (samps/sec).
%  cl.dt            Time domain bin width (ms).
%  cl.df            Frequency domain bin width (Hz).
%  cl.f_c95         95% confidence limit for Log spectral estimates.
%  cl.ch_c95        95% confidence limit for coherence.
%  cl.q_c95         95% confidence limit for cumulant density.
%  cl.a_c95         95% Confidence limits for impulse response (with s option).
%  cl.col_g         Column containing log gain (with s option).
%  cl.col_a         Column containing impulse response (with s option).
%  col_cova         Column containing Periodogram COV test on  input channel (with c option).
%  col_covb         Column containing Periodogram COV test on output channel (with c option).
%
%  sc column 1      f11.
%  sc column 2      f22.
%  sc column 3      f21 (Complex).
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%    Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Bloomfield, P. Fourier Analysis of Time Series: An Introduction.
%    2nd edition. Wiley, New York, 2000.
% 3. Nielsen J.B., Conway B.A., Halliday D.M., Perreault M-C. and Hultborn H.
%    Journal of Physiology, 569, 291-304, 2005.
%
% function [f,t,cl,sc] = sp2_fn2a(d1,d2,seg_samp_no,samp_rate,flags);
%
% NOTE: This routine is not intended to support analysis of raw data.
% It is intended as a support routine for the 2 channel spectral analysis functions:
% sp2_m1.m, sp2a_m1.m, sp2a2_m1.m. Refer to these functions for further details.

% PBMB     refers to above Progress in Biophysics article.
% JPHYSIOL refers to above  Journal of Physiology article.

[seg_size,seg_tot]=size(d1); % Calculate seg_size & seg_tot from data matrices.
samp_tot=sum(seg_samp_no);   % Total no of samples, Sum(T).

% Check if matrix of periodogram and cross periodogram coefficients required
flags.matrix=0;   % Default - No
if (flags.cov)    % Required in Periodogram COV test.
  flags.matrix=1;
end  

% Apply data window, if required.
if (flags.wind)
  switch (flags.wind)
    case 1       % 10%  split cosine taper applied  (5% at each end).
      p=0.1;
      u2=0.9377; % U2 factor (Bloomfield, P181).
    case 2       % 20%  split cosine taper applied (10% at each end).
      p=0.2;
      u2=0.8751;
    case 3       % 50%  split cosine taper applied (25% at each end).
      p=0.5;
      u2=0.6877;
    case 4       % 100%  full cosine taper applied (50% at each end).
      p=1;
      u2=0.3752;
  end

  % Calculate length of default taper.
  % Type 0, 1: data segment is length seg_size.
  % Type    2: data segment is length flags.seg_pts.
  if (flags.sp_type==2)
    tap_pts_tot=flags.seg_pts;
  else  
    tap_pts_tot=seg_size;
  end

  % Generate lower taper section for default taper.
  % No of points to taper at each end.
  tap_pts=round(tap_pts_tot*p/2);
  % Generate lower half of cosine over this number of points.
  taper(1:tap_pts,1)=0.5-0.5*cos(pi*(1:tap_pts)'/tap_pts);  % Bloomfield, (6.9).
  % Indices for upper data section to taper, in reverse order.
  tap_ind_up=tap_pts_tot:-1:tap_pts_tot-tap_pts+1;

  % Taper each segment
  for ind=1:seg_tot
    if (seg_samp_no(ind)==tap_pts_tot)
      % Use default taper.
      d1(1:tap_pts,ind) =d1(1:tap_pts,ind) .*taper; % Lower section.
      d1(tap_ind_up,ind)=d1(tap_ind_up,ind).*taper; % Upper section.
      d2(1:tap_pts,ind) =d2(1:tap_pts,ind) .*taper; % Lower section.
      d2(tap_ind_up,ind)=d2(tap_ind_up,ind).*taper; % Upper section.
      tap_pts_no(ind,1)=tap_pts;
    else
      % Have smaller number of data points - generate custom taper.
      % No of points to taper at each end.
      tap_pts_1=round(seg_samp_no(ind)*p/2);
      % Generate lower half of cosine over this number of points.
      taper_1=[];  % Clear any previous taper.
      taper_1(1:tap_pts_1,1)=0.5-0.5*cos(pi*(1:tap_pts_1)'/tap_pts_1); % Bloomfield, (6.9).
      % Indices for upper data section to taper in reverse order.
      tap_ind_up_1=seg_samp_no(ind):-1:seg_samp_no(ind)-tap_pts_1+1;
      d1(1:tap_pts_1,ind) =d1(1:tap_pts_1,ind) .*taper_1; % Lower section.
      d1(tap_ind_up_1,ind)=d1(tap_ind_up_1,ind).*taper_1; % Upper section.
      d2(1:tap_pts_1,ind) =d2(1:tap_pts_1,ind) .*taper_1; % Lower section.
      d2(tap_ind_up_1,ind)=d2(tap_ind_up_1,ind).*taper_1; % Upper section.
      tap_pts_no(ind,1)=tap_pts_1;
    end
  end

  % Normalisation factor for DFT to take account of tapering.
  u2_norm=sqrt(u2);            % Bloomfield, (9.11).
  fd1=fft(d1)/u2_norm;         % Take DFT across columns/segments ch 1, PBMB (4.1)/(4.2).
  fd2=fft(d2)/u2_norm;         % Take DFT across columns/segments ch 2, PBMB (4.1)/(4.2).

else
  % No tapering
  fd1=fft(d1);                 % Take DFT across columns/segments ch 1, PBMB (4.1)/(4.2).
  fd2=fft(d2);                 % Take DFT across columns/segments ch 2, PBMB (4.1)/(4.2). 
end
t_fac=2*pi*samp_tot;           % Normalization for weighted periodogram spectral estimates.
                               % JPHYSIOL (5), (6).

% Constuct 2D matrix of periodograms for individual segments.
if flags.matrix
  % Definition of indices to include - DC to f_n.
  per_ind=(1:seg_size/2+1)';
  for ind=1:seg_tot
    if flags.inv % Channel reference invert:  Channel 2 is input, Channel 1 is output.
      per_2d(:,1,ind)=abs(fd2(per_ind,ind).*fd2(per_ind,ind)/t_fac);   % Periodogram ch1
      per_2d(:,2,ind)=abs(fd1(per_ind,ind).*fd1(per_ind,ind)/t_fac);   % Periodogram ch2
      per_2d(:,3,ind)=fd1(per_ind,ind).*conj(fd2(per_ind,ind)/t_fac);  % Cross periodogram (complex)
    else         % Channels in default order: Channel 1 is input, Channel 2 is output.
      per_2d(:,1,ind)=abs(fd1(per_ind,ind).*fd1(per_ind,ind)/t_fac);   % Periodogram ch1
      per_2d(:,2,ind)=abs(fd2(per_ind,ind).*fd2(per_ind,ind)/t_fac);   % Periodogram ch2
      per_2d(:,3,ind)=fd2(per_ind,ind).*conj(fd1(per_ind,ind)/t_fac);  % Cross periodogram (complex)
    end  
  end
end

% Construct spectra based on average across segments
if flags.inv
  % Channel reference invert: Channel 2 is input, Channel 1 is output.
  f11=sum(abs(fd2.*fd2)/t_fac,2);   % Spectrum 1, PBMB (5.2), Mag squared for  input auto spectra.
  f22=sum(abs(fd1.*fd1)/t_fac,2);   % Spectrum 2, PBMB (5.2), Mag squared for output auto spectra.
  f21=sum(fd1.*conj(fd2)/t_fac,2);  % Cross spectrum (complex valued), PBMB (5.2).
else
  % Default order - Channel 1 is input, Channel 2 is output.
  f11=sum(abs(fd1.*fd1)/t_fac,2);   % Spectrum 1, PBMB (5.2), Mag squared for  input auto spectra.
  f22=sum(abs(fd2.*fd2)/t_fac,2);   % Spectrum 2, PBMB (5.2), Mag squared for output auto spectra.
  f21=sum(fd2.*conj(fd1)/t_fac,2);  % Cross spectrum (complex valued), PBMB (5.2).
end

deltaf=samp_rate/seg_size;        % Spacing of Fourier frequencies in Hz.

% Set line frequency in Hz.
line_freq=50;
% line_freq=60;  % Uncomment this line to set the line frequency to 60 Hz.

% Suppression of mains/line frequency - smooth out using adjacent values.
if flags.line      
  line_ind=round(line_freq/deltaf)+1;        % NB Index 1 is DC.
  f11(line_ind)=0.5*(f11(line_ind-2)+f11(line_ind+2));    % Spectrum ch 1.
  f11(line_ind-1)=0.5*(f11(line_ind-2)+f11(line_ind-3));
  f11(line_ind+1)=0.5*(f11(line_ind+2)+f11(line_ind+3));
  f22(line_ind)=0.5*(f22(line_ind-2)+f22(line_ind+2));    % Spectrum ch 2.
  f22(line_ind-1)=0.5*(f22(line_ind-2)+f22(line_ind-3));
  f22(line_ind+1)=0.5*(f22(line_ind+2)+f22(line_ind+3));
  f21(line_ind)=0.5*(f21(line_ind-2)+f21(line_ind+2));    % Cross spectrum.
  f21(line_ind-1)=0.5*(f21(line_ind-2)+f21(line_ind-3));
  f21(line_ind+1)=0.5*(f21(line_ind+2)+f21(line_ind+3));
% Smooth elements in upper hermetian section of cross spectral estimate.
% This data used in ifft() to generate cumulant. Data is complex conjugate.
  f21(seg_size-line_ind+2)=conj(f21(line_ind));
  f21(seg_size-line_ind+3)=conj(f21(line_ind-1));
  f21(seg_size-line_ind+1)=conj(f21(line_ind+1));
end

% Estimate cumulant density using inverse DFT of cross spectrum.
deltat=1000.0/samp_rate;    % dt in msec.

cov=ifft(f21);              % Inverse DFT.

% Construct output time domain matrix t.
% Column 1 of t matrix is time in msec. Range (-S/2)*dt to (S/2-1)*dt.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;

% Column 2 of t matrix is cumulant, shifted by S/2 so that time zero is in centre.
% 2pi/S factor is 2*pi, since ifft routine includes 1/S term.
t([seg_size/2+1:seg_size,1:seg_size/2],2)=real(cov(1:seg_size))*2*pi; % PBMB (5.9).

% Estimate variance of cumulant density estimate.
var_fac=4*pi*pi/(seg_size*samp_tot);                           % Factor (2pi/S)(2pi/R).
q_var=var_fac*2*sum(f11(1:seg_size/2+1).*f22(1:seg_size/2+1)); % PBMB (6.10).

% Optional impulse response calculation.
if flags.gain
  % Estimate impulse response using inverse DFT.
  f21_gain(1,1)=0;    % Disregard DC value, may well be zero.
  % Calculate gain, PBMB (11.2), f11 is input spectrum.
  f21_gain(2:seg_size,1)=f21(2:seg_size)./f11(2:seg_size);

  a=ifft(f21_gain);

  % Append impulse response to t matrix, shifted by T/2 so that time zero is in centre.
  % 1/T factor not included, since ifft routine includes 1/T term.
  col_a=size(t,2)+1;
  t([seg_size/2+1:seg_size,1:seg_size/2],col_a)=real(a(1:seg_size)); % PBMB (11.5).

  % Calculate variance of impulse response estimate.
  a_var_fac=1/(seg_size*samp_tot);                                 % Factor (1/T)(1/R).
  a_var=a_var_fac*2*sum(f22(2:seg_size/2+1)./f11(2:seg_size/2+1)); % PBMB (12.4).
end

% Reduce arrays to lower hermitian section only
f11=f11(1:seg_size/2+1);
f22=f22(1:seg_size/2+1);
f21=f21(1:seg_size/2+1);

% Apply smoothing using hanning filter.
var_smooth=1.0;     % Default - no smoothing correction.
if (flags.han==1)
  f11=han(f11);
  f22=han(f22);
  f21=han(f21);
  if flags.gain
    % Recalculate gain from smoothed spectra. PBMB (11.2), f11 is input spectrum.
    % Calculated only over lower hermitian section.
    f21_gain(2:seg_size/2+1,1)=f21(2:seg_size/2+1)./f11(2:seg_size/2+1);
  end
  var_smooth=0.375; % Correction for variance of smoothed estimates, Bloomfield (9.7).
end

% Construct output spectral matrix f.
f_index=(2:seg_size/2+1)';               % Indexing for output, DC component not output.
% Set minimum frequency output value for type two analysis, which takes account of zero padding.
if (flags.sp_type ==2)
  f_res=samp_rate/flags.seg_pts;         % Frequency resolution for type 2 analysis.
  % Calculate starting index for f matrix, with resolution deltaf.
  % Gives first frequency ordinate above f_res, +1 to allow for DC vlaue in first bin.
  f_index_start=ceil(f_res/deltaf)+1;
  f_index=(f_index_start:seg_size/2+1)'; % Indexing for type 2 analysis.
end
f(:,1)=(f_index-1)*deltaf;    % Column 1 - frequencies in Hz.
f(:,2)=log10(f11(f_index));   % Column 2 - Log spectrum ch 1.
f(:,3)=log10(f22(f_index));   % Column 3 - Log spectrum ch 2.
                              % Column 4 - Coherence, PBMB (5.5).
f(:,4)=abs(f21(f_index)).*abs(f21(f_index))./(f11(f_index).*f22(f_index));
f(:,5)=angle(f21(f_index));   % Column 5 - Phase, PBMB (5.7).

% Optional gain calculation.
if flags.gain
  % Append log magnitude of gain to f matrix, PBMB (11.3)
  col_g=size(f,2)+1;
  f(:,col_g)=log10(abs(f21_gain(f_index)));
end

% Optional periodogram COV estimates
if flags.cov
  col_cova=size(f,2)+1;
  col_covb=size(f,2)+2;
  f(:,col_cova)=std(per_2d(f_index,1,:),0,3)./mean(per_2d(f_index,1,:),3);  % COV of periodogram Ch 1
  f(:,col_covb)=std(per_2d(f_index,2,:),0,3)./mean(per_2d(f_index,2,:),3);  % COV of periodogram Ch 2
end

% Construct optional output spectral matrix sc, DC value is output.
if (nargout>3)
  sc(:,1)=f11(1:seg_size/2+1);  % Column 1 - f11
  sc(:,2)=f22(1:seg_size/2+1);  % Column 2 - f22
  sc(:,3)=f21(1:seg_size/2+1);  % Column 3 - f21 (Complex valued)
end

% Calculate effective number of segments, L', used to approximate variance of estimates.
seg_tot_var=1/sum((seg_samp_no(:,1)/samp_tot).*(seg_samp_no(:,1)/samp_tot)); % JPHYSIOL (8).
% Correct for smoothing.
seg_tot_var=seg_tot_var/var_smooth;

% Construct cl structure, confidence limits for parameter estimates.
cl.type=flags.sp_type;       % Analysis type.
cl.seg_size=seg_size;        % S.
cl.seg_tot=seg_tot;          % L.
cl.seg_tot_var=seg_tot_var;  % Effective no of segments (L').
cl.samp_tot=samp_tot;        % R.
cl.samp_rate=samp_rate;      % Sampling rate.
cl.dt=deltat;                % Delta t.
cl.df=deltaf;                % Delta f.
cl.f_c95=0.8512*sqrt(1/seg_tot_var);  % 95% Confidence limit for spectral estimates, PBMB (6.2).
% N.B. Confidence interval for log plot of spectra is TWICE this value.
cl.ch_c95=1-0.05^(1/(seg_tot_var-1)); % 95% Confidence limit for coherence, PBMB (6.6).
cl.q_c95=1.96*sqrt(q_var);            % 95% Confidence limits for cumulant, PBMB (6.11).

% Extra values for optional gain and impulse response estimates.
if flags.gain
  cl.a_c95=1.96*sqrt(a_var); % 95% Confidence limits for impulse response, PBMB (12.5).
  cl.col_g=col_g;            % Column in f matrix containing gain.
  cl.col_a=col_a;            % Column in t matrix containing impulse response.
end

if flags.cov
  cl.cov_c95=1.96*sqrt(3/(2*seg_tot));
  cl.col_cova=col_cova;
  cl.col_covb=col_covb;
end

return

%------------------------------------------------------------------------------
function [dat_out] = han(dat_in);
% function [dat_out] = han(dat_in);
%
% Function to smooth data vector using Hanning filter with coefficients: (1/4, 1/2, 1/4).
% End points smoothed using two point filter: (1/2, 1/2).

a=1;
b=[0.25;0.5;0.25];
dat_pts=length(dat_in);

% Smooth using function: filter.
dat_out=filter(b,a,dat_in);

% Shift one place to cancel out one sample delay in filter function.
dat_out=[dat_out(2:dat_pts);0];

% Smooth end points using two point filter: (1/2, 1/2).
dat_out(1)=0.5*(dat_in(1)+dat_in(2));
dat_out(dat_pts)=0.5*(dat_in(dat_pts)+dat_in(dat_pts-1));
return

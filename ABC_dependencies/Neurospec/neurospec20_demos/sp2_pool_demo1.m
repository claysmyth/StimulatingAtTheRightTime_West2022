% Script sp2_pool_demo1.m
% Script to demonstrate use of NeuroSpec2.0
%  Pooled spectra and pooled coherence analysis
%
% Copyright 2008, David M. Halliday.
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

%--------------------------------------------------------------------------

% Number of sets of data to pool.
pool_tot=10;

% Number of samples in each set.
samp_set_tot=11000;
% For purposes of analysis: assume sampling rate is 1000/sec
samp_rate=1000;

% Define filter parameters
% Filtering done using 11 point moving average window.
a=1;
b=(ones(1,11)/11)';

% Set power for segment length as 10.
% T=2^10 = 1024.
seg_pwr=10;

% No options used here
opt_str='';

% Loop round all data sets
for ind=1:pool_tot

  % Generate data set
  dat=randn(samp_set_tot,2);
  % Filter column 1.
   dat_filt=filter(b,a,dat(:,1));
  % Column 3 is average of filtered version of column 1 and column 2.
  dat(:,3)=0.5*(dat_filt+dat(:,2));

  % Last set - add in additional sine component: Magnitude 0.2, frequency 10 Hz.
  if (ind==pool_tot)
    freq=10;
    t=1/samp_rate*(1:samp_set_tot)';
    sin_dat=0.2*sin(2*pi*t*freq);
     dat(:,1)=dat(:,1)+sin_dat;
    dat(:,3)=dat(:,3)+sin_dat;
  end  

  % Normalise variance
  % Could also do this using 'n' option with sp2a2_m1
  dat(:,1)=dat(:,1)/std(dat(:,1));
  dat(:,3)=dat(:,3)/std(dat(:,3));

  % Process columns 1 and 3, these are Correlated & Filtered
  [f1(:,:,ind),t1(:,:,ind),cl1(ind),sc1(:,:,ind)] = sp2a2_m1(0,dat(:,1),dat(:,3),samp_rate,seg_pwr,opt_str);
  cl1(ind).what=['Set: ',num2str(ind)];

  % Pooled analysis
  if (ind==1)
    % Separate call for first set, creates new pooled analysis.
    [plf1,plv1]=pool_scf(sc1(:,:,ind),cl1(ind)); 
  else
    % Pass pooled variables as arguments for further sets.
    [plf1,plv1]=pool_scf(sc1(:,:,ind),cl1(ind),plf1,plv1);
  end  

end

% Plotting parameters
freq=75;
ch_max=1;
lag_tot=100;
lag_neg=50;
chi_max=0;% Will auto scale

% Process pooled spectral coefficients & plot
[f2,t2,cl2]=pool_scf_out(plf1,plv1);
figure
cl2.what='Pooled analysis, 10 sets';
psp2_pool6(f2,t2,cl2,freq,lag_tot,lag_neg,ch_max,chi_max)

figure
psp2_pool8(f2,t2,cl2,freq,lag_tot,lag_neg,ch_max,chi_max)

chi_max=125;  % All to same scale for comparison.
% Plots including three chi-squared tests (spectra, spectra, coherence).
% Eight plot version - No time domain
figure
psp2_pool8_chi3(f2,cl2,freq,ch_max,chi_max)

% Nine plot version - Includes time domain
figure
psp2_pool9_chi3(f2,t2,cl2,freq,lag_tot,lag_neg,ch_max,chi_max)


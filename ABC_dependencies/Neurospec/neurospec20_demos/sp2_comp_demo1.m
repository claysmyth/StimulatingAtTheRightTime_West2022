% Script sp2_comp_demo1.m
% Script to demonstrate use of NeuroSpec2.0
%  Comparison of spectra and comparison of coherence routines
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
% Generate data: Two channels with 10^5 points.
% Normally distributed, zero mean, unit variance: N(0,1).
dat=randn(100000,2);

% Generated filtered version of column 1.
% Filtering done using 11 point moving average window.
a=1;
b=(ones(1,11)/11)';
dat_filt=filter(b,a,dat(:,1));
% Column 3 is average of filtered version of column 1 and column 2.
dat(:,3)=0.5*(dat_filt+dat(:,2));

dat1=randn(100000,2);
% Generated filtered version of column 1.
% Filtering done using 21 point moving average window.
a=1;
b=(ones(1,21)/21)';
dat_filt=filter(b,a,dat1(:,1));
% Column 4 is average of filtered version of column 1 and column 2.
dat1(:,3)=0.5*(dat_filt+dat1(:,2));

% For purposes of analysis: assume sampling rate is 1000/sec
samp_rate=1000;

% Set power for segment length as 10.
% T=2^10 = 1024.
seg_pwr=10;

% No options used here
opt_str='';

% Process columns 1 and 3, these are Correlated & Filtered
[f1,t1,cl1,sc1]=sp2a2_m1(0,dat(:,1),dat(:,3),samp_rate,seg_pwr,opt_str);
cl1.what='Correlated & filtered (11pt)';
figure
freq=100;     % Parameters for plotting.
lag_tot=100;
lag_neg=50;
ch_max=1;
% Plot results.
psp2(f1,t1,cl1,freq,lag_tot,lag_neg,ch_max)

% Process columns 1 and 4, these are Correlated & Filtered
[f2,t2,cl2,sc2]=sp2a2_m1(0,dat1(:,1),dat1(:,3),samp_rate,seg_pwr,opt_str);
cl2.what='Correlated & filtered (21pt)';
figure
psp2(f2,t2,cl2,freq,lag_tot,lag_neg,ch_max)

% Comparison of spectra test.
[f3,cl3]=sp2_compf(sc1,cl1,2,sc2,cl2,2);
cl3.what='Log ratio test: f11/f22';
% Plot spectral estimates and log ratio test
freq=200;
figure
subplot(3,1,1)
psp_fb1(f1,cl1,freq)
subplot(3,1,2)
psp_fb1(f2,cl2,freq)
subplot(3,1,3)
psp_compf1(f3,cl3,freq)

% Comparison of coherence test
[f4,cl4]=sp2_compcoh(sc1,cl1,sc2,cl2);
cl4.what='Difference of coherence: Rab-Rcd';
% Plot coherence estimates and difference of coherence
freq=200;
ch_max=.8;
max_cmp=1;
figure
subplot(3,1,1)
psp_ch1(f1,cl1,freq,ch_max)
subplot(3,1,2)
psp_ch1(f2,cl2,freq,ch_max)
subplot(3,1,3)
psp_compcoh1(f4,cl4,freq,max_cmp)

% Script sp2_type2_demo1.m
% Script to demonstrate use of NeuroSpec 2.0
% TYPE 2 spectral analysis and plotting routines.
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
% Demonstration of TYPE 2 (time-frequency) analysis:
%  Common sine waves embedded into two independent N(0,1) noise signals.

%Set no of trials
trial_tot=50;

% Peak amplitude of sin waves. RMS magnitude is 0.707 times this.
sin10_mag=0.251; % This -15dB
sin25_mag=0.251; % This -15dB

% Each trial length is 1000 ms, 1000 points at 1 ms sampling interval.
trial_length=1000;
samp_rate=1000;
%No of samples in complete record.
samp_tot=trial_tot*trial_length;

% Generate 2 columns of N(0,1) data, @ 1ms sampling
dat=randn(samp_tot,2);

% Generate triggers: 1 per second over complete record.
% Triggers start at first sample in record.
trig=(1:trial_length:samp_tot)';

% Generate sin wave test data - each sin wave section has duration 100 ms.
% 100 ms time marker @ 1 ms
t=(0:0.001:0.099)';
% Generate 100 ms of 10 Hz sin wave, this is 1 cycle.
sin10=sin(2*pi*10*t);
% Generate 100 ms of 25 Hz sin wave, this is 2.5 cycles.
sin25=sin(2*pi*25*t);

% Offsets within each section for sin waves, specified as samples
sin10_offset=300;  % 10 Hz starts 300 samples into each section
sin25_offset=500;  % 25 Hz starts 500 samples into each section

% Add in sin waves.
for ind=1:length(trig)

  trial_start=trig(ind);
  % Add in 10 Hz component
  sin10_start=trial_start+sin10_offset;
  sin10_stop=sin10_start+length(sin10)-1;
  dat(sin10_start:sin10_stop,1)=dat(sin10_start:sin10_stop,1)+sin10_mag*sin10;
  dat(sin10_start:sin10_stop,2)=dat(sin10_start:sin10_stop,2)+sin10_mag*sin10;

  % Add in 25 Hz component
  sin25_start=trial_start+sin25_offset;
  sin25_stop=sin25_start+length(sin25)-1;
  dat(sin25_start:sin25_stop,1)=dat(sin25_start:sin25_stop,1)+sin25_mag*sin25;
  dat(sin25_start:sin25_stop,2)=dat(sin25_start:sin25_stop,2)+sin25_mag*sin25;
end  

%------------------------------------------------------------------------------
% Analysis

% Offset values (in samples) for TYPE 2 time-frequency analysis.
offset=(0:50:750)';

% No of data points in each segment (samples)
seg_pts=250;
% DFT length for time-freq analysis
seg_pwr=8;

% Options
opt_str='';

[f,t,cl]=sp2a2_m1(2,dat(:,1),dat(:,2),trig,offset,seg_pts,samp_rate,seg_pwr,opt_str);
cl(1).what='N(0,1) + 10, 25 Hz';

% Plot using new routine, using default no of contours, and default line flag.
% psp2_tf(f,t,cl,freq,lag_tot,lag_neg,t_inc,n_contour,line_flag,coh_contour)
freq=50;
lag_tot=150;
lag_neg=75;
psp2_tf(f,t,cl,freq,lag_tot,lag_neg)

% Plot two sections through time frequency plane at fixed offset
sect_1=7;   % Index of 10 Hz section
sect_2=11;  % Index of 25 Hz section

freq=100;
lag_tot=200;
lag_neg=100;
ch_max=0.3;

figure
psp2(f(:,:,sect_1),t(:,:,sect_1),cl(sect_1),freq,lag_tot,lag_neg,ch_max,[cl(1).what,' Offset: ',num2str(cl(sect_1).offset)])
figure
psp2(f(:,:,sect_2),t(:,:,sect_2),cl(sect_2),freq,lag_tot,lag_neg,ch_max,[cl(1).what,' Offset: ',num2str(cl(sect_2).offset)])

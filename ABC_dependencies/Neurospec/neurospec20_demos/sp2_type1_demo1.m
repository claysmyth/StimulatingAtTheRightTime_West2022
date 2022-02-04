% Script sp2_type1_demo1.m
% Script to demonstrate use of NeuroSpec 2.0
% TYPE 1 spectral analysis and plotting routines.
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
% Demonstration of TYPE 1 analysis:
%  Analysis using distinct sections from single record.


%--------------------------------------------------------------------------
% Section 1: Usage of sp2a2_m1.m  with Type 1 analysis

% Generate some data: Two channels with 10^5 points.
% Normally distributed, zero mean, unit variance: N(0,1).
dat=randn(100000,2);
% Add third column which is average of first two.
dat(:,3)=0.5*(dat(:,1)+dat(:,2));

% For purposes of analysis: assume sampling rate is 1000/sec
samp_rate=1000;

% Set power for segment length as 10.
% Thus T=2^10 = 1024.
seg_pwr=10;
seg_size=2^seg_pwr;

% Generate trig_times and duration for type 1 analysis, with random durations.
samp_tot=length(dat);   % No of samples available for analysis
trig_count=0;           % Counts triggers
samp_count=0;           % Counts samples
samp_min=50;            % Minimum number of samples in each segment
trig_times=[];          %  Vector for start sample of each segment
duration=[];            % Vector for duration of each segment
while samp_count<samp_tot
  seg_pts=samp_min+round(rand*(seg_size-samp_min)); % Try this number of points in segment
  if (samp_count+seg_pts>samp_tot)                  % If too many - shorten
    seg_pts=samp_tot-samp_count;                    % Reduce to fit length of data
  end  
  trig_count=trig_count+1;                          % New trigger
  trig_times(trig_count,1)=samp_count+1;            % Trigger start is 1 past previous point
  duration(trig_count,1)=seg_pts;                   % Set duration
  samp_count=samp_count+seg_pts;                    % Update sample count
end
% Now have valid set of triggers and durations for Type 1 analysis.
disp(['Type 1 - No of triggers: ',num2str(length(trig_times))])

% No options used here
opt_str='';

% Process columns 1 and 2, these are Uncorrelated
[f,t,cl]=sp2a2_m1(1,dat(:,1),dat(:,2),trig_times,duration,samp_rate,seg_pwr,opt_str);
cl.what='N(0,1) uncorrelated data.';
figure
% The sampling rate defines the nyquist frequency as 500 Hz.
freq=500;     % Parameters for plotting
lag_tot=500;
lag_neg=250;
ch_max=0.5;
psp2(f,t,cl,freq,lag_tot,lag_neg,ch_max)

% Process columns 1 and 3, these are Correlated (column 1 contributed to column 3)
[f1,t1,cl1]=sp2a2_m1(1,dat(:,1),dat(:,3),trig_times,duration,samp_rate,seg_pwr,opt_str);
cl1.what='Correlated data';
figure
freq=500;     % Parameters for plotting
lag_tot=50;
lag_neg=25;
ch_max=1;
psp2(f1,t1,cl1,freq,lag_tot,lag_neg,ch_max)

% Column 4 is delayed version of column3. Delay of 10 samples.
dat(1:10,4)=randn(10,1);         % First 10 samples generated independently
dat(11:100000,4)=dat(1:99990,3); % Subsequent samples in column 4 copied from column 3.

% Process columns 2 and 4, these are Correlated & Delayed
%  Column 2 contributed to column 3, column 4 is delayed version of column 3.
[f2,t2,cl2]=sp2a2_m1(1,dat(:,2),dat(:,4),trig_times,duration,samp_rate,seg_pwr,opt_str);
cl2.what='Correlated & delayed data 1';
figure
psp2(f2,t2,cl2,freq,lag_tot,lag_neg,ch_max)

% Generated filtered version of column 1.
% Filtering done using 11 point moving average window.
a=1;
b=(ones(1,11)/11)';
dat_filt=filter(b,a,dat(:,1));

% Column 5 is average of filtered version of column 1 and column 2.
dat(:,5)=0.5*(dat_filt+dat(:,2));

% Process columns 1 and 5, these are Correlated & Filtered
% Include System Identification analysis in options string
opt_str='s';

[f3,t3,cl3]=sp2a2_m1(1,dat(:,1),dat(:,5),trig_times,duration,samp_rate,seg_pwr,opt_str);
cl3.what='Correlated & filtered data 1';
figure
freq=100;     % Reduced frequency range for filtered data
% Plot with psp2s, includes optional system identication parameter estimates:
%  Gain and Impulse response.
psp2s(f3,t3,cl3,freq,lag_tot,lag_neg,ch_max)

%--------------------------------------------------------------------------
% Section 2: Usage of sp2_m1.m 

% Generate some spike train data.
% Clear these each time.
clear int1 int2 sp sp1 sp2

% First spike train, Mean Interval=100 sampling bins, STD=10 sampling bins.
% Generate 1200 inter spike intervals.
int1=round(randn(1200,1)*10+100);
% Sum these to get spike/event times.
for ind=1:length(int1)
  sp(ind,1)=sum(int1(1:ind));
end
% Now copy to sp1 those spike times less than 100000
% This represents 100 seconds @ 1 ms sampling, same as time series data.
sp1(1:length(find(sp<100000)),1)=sp(1:length(find(sp<100000)));

% Second spike train, Mean Interval=100 sampling bins, STD=20 sampling bins.
% Generate 1200 inter spike intervals.
int2=round(randn(1200,1)*20+100);
% Sum these to get spike/event times.
for ind=1:length(int2)
  sp(ind,1)=sum(int2(1:ind));
end 
% Now copy to sp2 those spike times less than 100000
% This represents 100 seconds @ 1 ms sampling, same as time series data.
sp2(1:length(find(sp<100000)),1)=sp(1:length(find(sp<100000)));

% Now now have two independent spike trains both with a mean rate of 10 spikes/sec,
% the first with an STD of 10 sampling bins, the second with an STD of 20 sampling bins.

% Process these using sp2_m1
sec_tot=100;   % sp2_m1 requires no of seconds as an argument.
[f4,t4,cl4]=sp2_m1(1,sp1,sp2,trig_times,duration,sec_tot,samp_rate,seg_pwr,opt_str);
cl4.what='Independent Spike Trains';
figure
freq=75;     % Parameters for plotting
lag_tot=100;
lag_neg=50;
ch_max=0.5;
psp2(f4,t4,cl4,freq,lag_tot,lag_neg,ch_max)

%--------------------------------------------------------------------------
% Section 3: Usage of sp2a_m1.m

% Now have spike train (point process) and time series data.Do hybrid analysis using
%  spike train data from section 2 and time series data from section 1.

[f5,t5,cl5]=sp2a_m1(1,sp1,dat(:,1),trig_times,duration,samp_rate,seg_pwr,opt_str);
cl5.what='Hybrid analysis';
figure
psp2(f5,t5,cl5,freq,lag_tot,lag_neg,ch_max)

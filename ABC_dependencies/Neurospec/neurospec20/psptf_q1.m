function psptf_q1(t,cl,lag_tot,lag_neg,t_inc,n_contour,line_flag)
% function psptf_q1(t,cl,lag_tot,lag_neg,t_inc,n_contour,line_flag)
% Function to display time dependent cumulant estimate in current subplot window.
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
% Inputs
%    t              Time matrix from type 2 analysis, over range of offset values.
%    cl             cl structure from type 2 analysis.
%    lag_tot        Total lag range for time domain including -ve lags (ms).
%    lag_neg        Negative lag range for time domain (ms).
%    t_inc          Optional increment for labels on time axis (0: auto scales).
%    n_contour      Optional number of contours to use (default 10).
%    colorbar_flag  Optional flag to control drawing of colorbar      - 0:No(default), 1:Yes.
%    line_flag      Optional flag to control drawing of contour lines - 0:No(default), 1:Yes.
%
% function psptf_q1(t,cl,lag_tot,lag_neg,t_inc,n_contour,line_flag)

% Check numbers of arguments. 
if (nargin<4)
  error('Not enough input arguments');
end
if (cl(1).type~=2)
  error('Not type 2 analysis')
end

% Defaults
n_seg=length(t(1,1,:));  % No of time slice segments.
if (n_seg<3)
  error('Less than 3 time offsets')
end
if (nargin<5)
  t_inc=0;
end
if (t_inc==0)
  if (n_seg<=10)
    t_inc=1;
    elseif (n_seg<=20)
      t_inc=2;
    elseif (n_seg<=50)
      t_inc=5;
    elseif (n_seg<=100)
      t_inc=10;
    else
      t_inc=25;
  end
end
if (nargin<6)
  n_contour=10;
end  
if (nargin<7)
  line_flag=0;
end  

bin_start=cl(1).seg_size/2+1-round(lag_neg/cl(1).dt);
bin_stop=bin_start+round(lag_tot/cl(1).dt)-1;
%Check lag range
[x,y]=size(t(:,:,1));
if (bin_start<1 | bin_start>=bin_stop | bin_stop>x)
  error('Error in requested lag range.');
end
bin_pts=bin_stop-bin_start+1;

% Set up x labels in   XTickLabel format.
dt=cl(1).dt;             % Sampling interval (ms). 
t_inc=round(t_inc);      % t_inc should be integer.
t_array=(1:t_inc:n_seg); % Array of indices for time ticks.
t_str=[];
for ind=1:t_inc:n_seg
  if (ind<n_seg)
    t_str=strcat(t_str,num2str(cl(ind).offset*dt),'|');
  else
    t_str=strcat(t_str,num2str(cl(ind).offset*dt));
  end  
end

%bin_zero=cl(1).seg_size/2+1;    % Index of zero lag

% Define values for cumulant lag ticks (ms)
if (lag_tot<=50)
  y_start=fix(-lag_neg/5);
  y_lag=(y_start*5:5:lag_tot-lag_neg);
  elseif (lag_tot<=100)
    y_start=fix(-lag_neg/10);
    y_lag=(y_start*10:10:lag_tot-lag_neg);
  elseif (lag_tot<=250)
    y_start=fix(-lag_neg/20);
    y_lag=(y_start*20:20:lag_tot-lag_neg);
  else
    y_start=fix(-lag_neg/50);
    y_lag=(y_start*50:50:lag_tot-lag_neg);
end
% Set up labels in   YTickLabel format.
y_str=[];
for ind=1:length(y_lag)
  if (ind<length(y_lag))
    y_str=strcat(y_str,num2str(y_lag(ind)),'|');
  else
    y_str=strcat(y_str,num2str(y_lag(ind)));
  end
end

% Calculate index in t array of y_lag values
y_val=(cl(1).seg_size/2+1+y_lag/dt);

% Reduce array for plotting.
tq(1:bin_pts,1:n_seg)=t(bin_start:bin_stop,2,1:n_seg);
% Correct y_val for reduced array.
y_val=y_val-bin_start+1;

[C,h]=contourf(tq,n_contour);
H1=gca;
if (line_flag==0)
  for ind=1:length(h)
    set(h(ind),'LineStyle','none')
  end
end
set (gca,'XTick',t_array)
set(gca,'XTickLabel',t_str)
set (gca,'YTick',y_val)
set(gca,'YTickLabel',y_str)
set (gca,'TickDir','out')
xlabel ('Offset (ms)')
ylabel ('Lag (ms)')
title([cl(1).what, ' - Time dependent cumulant'])

H=colorbar;
v=get(H);
axes(H)
for q_ind=1:n_seg
  q_c95u_vector=[v.XLim;cl(q_ind).q_c95,cl(q_ind).q_c95];
  line(q_c95u_vector(1,:),q_c95u_vector(2,:),'Color','k')
  q_c95l_vector=[v.XLim;-cl(q_ind).q_c95,-cl(q_ind).q_c95];
  line(q_c95l_vector(1,:),q_c95l_vector(2,:),'Color','k')
end
axes(H1)

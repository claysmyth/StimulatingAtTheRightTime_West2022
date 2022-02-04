function psptf_fb1(f,cl,freq,t_inc,n_contour,colorbar_flag,line_flag)
% function psptf_fb1(f,cl,freq,t_inc,n_contour,colorbar_flag,line_flag)
% Function to display time dependent output auto spectral estimate in current subplot window.
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
%    f              Frequency matrix from type 2 analysis, over range of offset values.
%    cl             cl structure from type 2 analysis.
%    freq           Frequency limit for plotting.
%    t_inc          Optional increment for labels on time axis (0: auto scales).
%    n_contour      Optional number of contours to use (default 10).
%    colorbar_flag  Optional flag to control drawing of colorbar      - 0:No(default), 1:Yes.
%    line_flag      Optional flag to control drawing of contour lines - 0:No(default), 1:Yes.
%
% function psptf_fb1(f,cl,freq,t_inc,n_contour,colorbar_flag,line_flag)

% Check numbers of arguments. 
if (nargin<3)
  error('Not enough input arguments');
end
if (cl(1).type~=2)
  error('Not type 2 analysis')
end

% Defaults
n_seg=length(f(1,1,:));  % No of time slice segments.
if (n_seg<3)
  error('Less than 3 time offsets')
end
if (nargin<4)
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
if (nargin<5)
  n_contour=10;
end  
if (nargin<6)
  colorbar_flag=0;
end  
if (nargin<7)
  line_flag=0;
end  

f_start=f(1,1,1);  % Starting frequency for plots.
df=cl(1).df;       % Frequency spacing of estimates.
% Define values for frequency ticks (Hz), dependent on frequency range.
if (freq<=50)
  f_val=(5:5:freq);
  elseif (freq<=100)
    f_val=(10:10:freq);
  elseif (freq<=250)
    f_val=(20:20:freq);
  else
    f_val=(50:50:freq);
end
% Set up labels in   YTickLabel format.
f_str=[];
for ind=1:length(f_val)
  if (ind<length(f_val))
    f_str=strcat(f_str,num2str(f_val(ind)),'|');
  else
    f_str=strcat(f_str,num2str(f_val(ind)));
  end
end
% Calculate position of indices that match frequencies in f_val.
f_array=(f_val-f_start)/df+1;
% No of points to plot in frequency axis.
freq_pts=round((freq-f_start)/cl(1).df);

dt=cl(1).dt;             % Sampling interval (ms). 
t_inc=round(t_inc);      % t_inc should be integer.
t_array=(1:t_inc:n_seg); % Array of indices for time ticks.
% Set up labels in   XTickLabel format.
t_str=[];
for ind=1:t_inc:n_seg
  if (ind<n_seg)
    t_str=strcat(t_str,num2str(cl(ind).offset*dt),'|');
  else
    t_str=strcat(t_str,num2str(cl(ind).offset*dt));
  end  
end

% Plot Time dependent log spectra of channel 2.
[C,h]=contourf(squeeze(f(1:freq_pts,3,1:n_seg)),n_contour);
if (line_flag==0)
  for ind=1:length(h)
    set(h(ind),'LineStyle','none')
  end
end  
set (gca,'XTick',t_array)
set(gca,'XTickLabel',t_str)
set (gca,'YTick',f_array)
set(gca,'YTickLabel',f_str)
set (gca,'TickDir','out')
xlabel('Offset (ms)')
ylabel ('Frequency (Hz)')
title([cl(1).what, ' Time dependent Spectra'])
if colorbar_flag
  colorbar
end  

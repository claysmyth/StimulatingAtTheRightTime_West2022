function psp_ch1cl(f,cl,freq,ch_max,plot_type,label)
% psp_ch1cl(f,cl,freq,ch_max,plot_type,label)
% function to plot coherence in current figure/subplot window
% includes pointwise confidence limits.
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
% f,cl      Output from spectral analysis routine.
% freq      Frequency limit for plotting (Hz).
% ch_max    Maximum value of y axis (Optional).
% plot_type Type of display - 0: Line width (default); 1: Line style; 2: Error bar.
% label     Optional title instead of cl.what.
%
%  psp_ch1cl(f,cl,freq,ch_max,plot_type,label)

if (nargin<3)
  error('Not enough input arguments')
end  

%freq_pts=round(freq/cl.df);
%Updated to support type 2 analysis
f_start=f(1,1,1);
freq_pts=round((freq-f_start)/cl.df)+1;

%Check freq range
[x,y]=size(f);
if (freq_pts>x)
  error('Requested frequency range too large.');
end
if (freq_pts<1)
  error('Requested frequency range too narrow.');
end

f_max=f(freq_pts,1);
% Estimate of upper 95% confidence limit for independence
ch_cl=[f_start,cl.ch_c95;f_max,cl.ch_c95];

% Upper and lower 95% pointwise confidence limits, based on Halliday et al. (1995), (6.5).
ch_c95up=tanh(atanh(sqrt(f(1:freq_pts,4)))+(1.96/sqrt(2*cl.seg_tot_var))).^2;
ch_c95lp=tanh(atanh(sqrt(f(1:freq_pts,4)))-(1.96/sqrt(2*cl.seg_tot_var))).^2;

% Convert absolute value of upper and lower 95% pointwise limits to relative distance
% for use with errorbar plotting.
ch_c95up_errbar=f(1:freq_pts,4)-ch_c95up(1:freq_pts,1);
ch_c95lp_errbar=f(1:freq_pts,4)-ch_c95lp(1:freq_pts,1);

if (nargin<5)
  plot_type=0;
end

switch plot_type
  case 0
  % Plot coherence estimate plus confidence limit for independence
    plot(f(1:freq_pts,1),f(1:freq_pts,4),'k-','LineWidth',2)
    hold on
    plot(ch_cl(:,1),ch_cl(:,2),'k--');
  % Point wise limits as thinner lines
    plot(f(1:freq_pts,1),ch_c95up(1:freq_pts,1),'k-');
    plot(f(1:freq_pts,1),ch_c95lp(1:freq_pts,1),'k-');
  case 1
  % Plot coherence estimate plus confidence limit for independence
    plot(f(1:freq_pts,1),f(1:freq_pts,4),'k-')
    hold on
    plot(ch_cl(:,1),ch_cl(:,2),'k--');
  % Point wise limits as dotted lines
    plot(f(1:freq_pts,1),ch_c95up(1:freq_pts,1),'k:');
    plot(f(1:freq_pts,1),ch_c95lp(1:freq_pts,1),'k:');
  case 2
  % Using error bar command
    errorbar(f(1:freq_pts,1),f(1:freq_pts,4),ch_c95lp_errbar(1:freq_pts,1),ch_c95up_errbar(1:freq_pts,1),'k')
    hold on
  % Fixed confidence limit for independence
    plot(ch_cl(:,1),ch_cl(:,2),'k--');
  otherwise
    error ('Unknown plot_type')
end
hold off

if (nargin<4)
  axis([0,freq,0,Inf]);
else
  axis([0,freq,0,ch_max]);
end  
xlabel('Frequency (Hz)')
if (nargin>5)
  title(label);
else
 title(['coh: ',cl.what]);
end  
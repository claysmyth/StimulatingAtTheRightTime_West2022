function psp_ph1cl(f,cl,freq,plot_type,label)
% psp_ph1cl(f,cl,freq,plot_type,label)
% function to plot phase in current figure/subplot window
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
%     Contact:  contact@neurospec.org
%
%  psp_ph1cl(f,cl,freq,plot_type,label)
%
% f,cl      Output from spectral analysis routine.
% freq      Frequency limit for plotting (Hz).
% plot_type Type of display - 0: Line width (default); 1: Line style; 2: Error bar.
% label     Optional title instead of cl.what.

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
% Co-ordinates of zero phase line
ph_zero=[f_start,0.0;f_max,0.0];

% Upper and lower 95% pointwise confidence limits, based on Halliday et al. (1995), (6.8).
ph_c95up=f(:,5)+1.96*sqrt((1/(2*cl.seg_tot_var))*(1./f(:,4)-1));
ph_c95lp=f(:,5)-1.96*sqrt((1/(2*cl.seg_tot_var))*(1./f(:,4)-1));

% Convert absolute value of upper and lower 95% pointwise limits to relative distance
% for use with errorbar plotting.
ph_c95up_errbar=f(1:freq_pts,5)-ph_c95up(1:freq_pts,1);
ph_c95lp_errbar=f(1:freq_pts,5)-ph_c95lp(1:freq_pts,1);

if (nargin<4)
  plot_type=0;
end

switch plot_type
  case 0
  % Plot phase
    plot(f(1:freq_pts,1),f(1:freq_pts,5),'k-','LineWidth',2);
    hold on
  % Point wise limits as thinner lines
    plot(f(1:freq_pts,1),ph_c95up(1:freq_pts,1),'k-');
    plot(f(1:freq_pts,1),ph_c95lp(1:freq_pts,1),'k-');
    plot(ph_zero(:,1),ph_zero(:,2),'k--');
  case 1
  % Plot phase
    plot(f(1:freq_pts,1),f(1:freq_pts,5),'k-');
    hold on
  % Point wise limits as dotted lines
    plot(f(1:freq_pts,1),ph_c95up(1:freq_pts,1),'k:');
    plot(f(1:freq_pts,1),ph_c95lp(1:freq_pts,1),'k:');
    plot(ph_zero(:,1),ph_zero(:,2),'k--');
  case 2
  % Using error bar command
    errorbar(f(1:freq_pts,1),f(1:freq_pts,5),ph_c95lp_errbar(1:freq_pts,1),ph_c95up_errbar(1:freq_pts,1),'k')
    hold on
    plot(ph_zero(:,1),ph_zero(:,2),'k--');
  otherwise
    error ('Unknown plot_type')
end
hold off

axis([0,freq,-pi,pi]);
xlabel('Frequency (Hz)')
if (nargin>4)
  title(label);
else
 title(['ph: ',cl.what]);
end  
function psp2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)
% psp2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)
% Function to plot output of 2 channel spectral analysis,
%  from type 0, type 1 or single offset type 2.
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
% f,t,cl   Output from spectral analysis routine.
% freq     Frequency limit for plotting (Hz).
% lag_tot  Total lag range for time domain including -ve lags (ms).
% lag_neg  Negative lag range for time domain (ms).    
% ch_max   Upper limit for coherence estimte.
% label    Label for subplots, used in preference to cl.what
%
% psp2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)

% Check numbers of arguments. 
if (nargin<7)
  error('Not enough input arguments');
end

psp_fa1(f,cl,freq,label)

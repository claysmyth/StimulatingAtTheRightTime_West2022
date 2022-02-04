% NeuroSpec
% Version 2.0  29-Feb-2008
%
%Core routines
%  sp2_m1          - Analysis of two spike trains.
%  sp2a_m1         - Analysis of one spike train and one time series.
%  sp2a2_m1        - Analysis of two time series.
%  sp2_fn2a        - Function called by routines (Two channel weighted periodogram analysis).
%
%Comparison of spectra and comparison of coherence
%  sp2_compf       - Comparison of spectra test.
%  sp2_compcoh     - Comparison of coherence test.
%
%Pooled spectra and coherence
%  pool_scf        - Create or update pooled spectral coe±cients.
%  pool_scf_out    - Convert pooled spectral coe±cients to form suitable for plotting.
%
%Main Plotting Routines
%  psp2            - Plots Type 0, Type 1 or Type 2 (single offset) analysis from core routines.
%  psp2s           - As psp2, includes optional system identification parameters.
%  psp2_cov        - As psp2, includes optional Periodogram COV test for each channel.
%  psp2_tf         - Plotting of Type 2 time-frequency analysis over range of offset values.
%  psp_compf1      - Plots results of comparison of spectra test.
%  psp_compcoh1    - Plots results of comparison of coherence test.
%  psp2_pool6      - Plots results of pooled analysis.
%  psp2_pool8      - Plots results of pooled analysis (includes 2 additional panels).
%  psp2_pool8_chi3 - Plots results of pooled analysis, includes 3 chi-squared tests (no time domain).
%  psp2_pool9_chi3 - Plots results of pooled analysis, includes 3 chi-squared tests (with time domain).
%
%Other Plotting Routines (used in)
%  psp_fa1         - Plots input log spectral estimate (psp2, psp2s, psp2 pool*)
%  psp_fb1         - Plots output log spectral estimate (psp2, psp2s, psp2 pool*)
%  psp_fcova1      - Plots Periodogram COV test for input spectral estimate.
%  psp_fcovb1      - Plots Periodogram COV test for output spectral estimate.
%  psp_ch1         - Plots coherence estimate (psp2, psp2s, psp2 pool*)
%  psp_ph1         - Plots phase estimate (psp2, psp2s, psp2 pool*)
%  psp_q1          - Plots cumulant density estimate (psp2, psp2s, psp2 pool*)
%  psp_g1          - Plots log gain estimate (psp2s)
%  psp_a1          - Plots impulse response estimate (psp2s)
%  psp_chi1        - Plots difference of coherence test (psp2 pool*)
%  psp_fchia1      - Plots difference of input spectra test (psp2 pool8 chi3, psp2 pool9 chi3)
%  psp_fchib1      - Plots difference of output spectra test (psp2 pool8 chi3, psp2 pool9 chi3)
%  psp_chplf1      - Plots alternative pooled coherence estimate (psp2 pool8)
%  psp_chst1       - Plots histogram of significant coherence values in population (psp2 pool*)
%  psptf_fa1       - Plots time dependent input log spectral estimate (psp2 tf)
%  psptf_fb1       - Plots time dependent output log spectral estimate (psp2 tf)
%  psptf_ch1       - Plots time dependent coherence estimate (psp2 tf)
%  psptf_ph1       - Plots time dependent phase estimate (psp2 tf)
%  psptf_q1        - Plots time dependent cumulant density estimate (psp2 tf)
%
%Plotting of individual estimates with point-wise confidence limits
%  psp_ch1cl       - Plots coherence estimate.
%  psp_ph1cl       - Plots phase estimate.
%  psp_g1cl        - Plots log gain estimate.
%
%Demonstration scripts
%  sp2_type0_demo1 - Demonstrates Type 0 analysis.
%  sp2_type1_demo1 - Demonstrates Type 1 analysis.
%  sp2_type2_demo1 - Demonstrates Type 2 analysis.
%  sp2_comp_demo1  - Demonstrates comparison of spectra and comparison of coherence.
%  sp2_pool_demo1  - Demonstrates pooled spectral and coherence analysis.

%Copyright (C) 2008, David M. Halliday.
%This file is part of NeuroSpec.
%
%   NeuroSpec is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   NeuroSpec is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with NeuroSpec; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%   NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org

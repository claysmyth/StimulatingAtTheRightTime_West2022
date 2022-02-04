function [f,t,cl,sc] = sp2_m1(sp_type,sp1,sp2,varargin);
% Type 0: [f,t,cl,sc] =  sp2_m1(0,sp1,sp2,sec_tot,samp_rate,seg_pwr,opt_str)
% Type 1: [f,t,cl,sc] =  sp2_m1(1,sp1,sp2,trig_times,duration,sec_tot,samp_rate,seg_pwr,opt_str)
% Type 2: [f,t,cl,sc] =  sp2_m1(2,sp1,sp2,trig_times,offset,seg_pts,sec_tot,samp_rate,seg_pwr,opt_str)
% Function to calculate spectra, coherence, phase & cumulant for 2 spike trains.
% Type 0, 1, 2 analysis - using weighted periodogram estimates.
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
% Input arguments
%  sp_type    Analysis type (0, 1, 2).
%  sp1        Channel 1   (input) spike times - specify as integer no of sampling intervals,
%                                               ordered in ascending order.
%  sp2        Channel 2  (output) spike times - specify as integer no of sampling intervals,
%                                               ordered in ascending order.
%
%  Additional arguments
%  Type 0:
%   sec_tot    Length of data set in seconds.
%   samp_rate  Sampling rate (samples/sec).
%   seg_pwr    Segment length - specified as power of 2.
%   opt_str    Options string.
%
%  Type 1:
%   trig_times Vector: List of trigger times defining start of each data segment (in samples).
%   duration   Vector: Duration of each data segment (in samples).
%   sec_tot    Length of data set in seconds.
%   samp_rate  Sampling rate (samples/sec).
%   seg_pwr    Segment length - specified as power of 2.
%   opt_str    Options string.
%
%  Type 2:
%   trig_times Vector: List of trigger times defining start of each data segment (in samples).
%   offset     Vector: List of offset values from trigger times to start of each data segment
%                      A separate analysis is done for each offset value.
%   seg_pts    Scalar: fixed number of data points in each segment (in samples).
%   sec_tot    Length of data set in seconds.
%   samp_rate  Sampling rate (samples/sec).
%   seg_pwr    Segment length - specified as power of 2.
%   opt_str    Options string.
%
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%  sc   optional matrix for spectral coefficients.
%
% Input Options
%  c  Include Periodogram COV test (Not valid for type 1 analysis).
%  h  Additional smoothing of spectra using hanning filter.
%  i  Invert channel reference for phase and cumulant density.
%  m  Mains/line frequency suppression.
%  s  System Identification: Estimates gain and impulse response.
%  w  Use data window  - this options requires an argument, valid arguments are:[1, 2, 3, 4].
%                      w1: Split cosine taper applied to each segment,  5% tapered at each end.
%                      w2: Split cosine taper applied to each segment, 10% tapered at each end.
%                      w3: Split cosine taper applied to each segment, 25% tapered at each end.
%                      w4: Full  cosine taper applied to each segment, 50% tapered at each end.
%
% Options examples:
%  to set all options on (5% split taper), set opt_str='c h i m s w1'
%
% Output parameters
%  f column 1       frequency in Hz.
%  f column 2       Log input/sp1  spectrum.
%  f column 3       Log output/sp2 spectrum.
%  f column 4       Coherence.
%  f column 5       Phase.
%  f column 6       Log of gain magnitude (with s option).
%  f column 6 or 7  Periodogram COV test  input channel (with c option).
%  f column 7 or 8  Periodogram COV test output channel (with c option).
%
%  t column 1       Lag in ms.
%  t column 2       Cumulant density.
%  t column 3       Impulse response (with s option).
%
%  cl.type          Analysis type (0, 1, 2)
%  cl.seg_size      Segment length.
%  cl.seg_tot       Number of segments.
%  cl.seg_tot_var   Effective no of segments, used to calculate confidence limits.
%  cl.samp_tot      Number of samples analysed.
%  cl.samp_rate     Sampling rate of data (samps/sec).
%  cl.dt            Time domain bin width (ms).
%  cl.df            Frequency domain bin width (Hz).
%  cl.f_c95         95% confidence limit for Log spectral estimates.
%  cl.ch_c95        95% confidence limit for coherence.
%  cl.q_c95         95% confidence limit for cumulant density.
%  cl.a_c95         95% confidence limit for impulse response      (with s option).
%  cl.col_g         column in f matrix containing gain estimate    (with s option).
%  cl.col_a         column in t matrix containing impulse response (with s option).
%  col_cova         Column containing Periodogram COV test on  input channel (with c option).
%  col_covb         Column containing Periodogram COV test on output channel (with c option).
%  cl.seg_pts       Data points in segment       (Type 2 only).
%  cl.offset        Offset value(s) for analysis (Type 2 only).
%  cl.N1            N1, No of events/spikes in sp1. 
%  cl.N2            N2, No of events/spikes in sp2.
%  cl.P1            P1, mean intensity sp1.
%  cl.P2            P2, mean intensity sp2.
%  cl.opt_str       Copy of options string.
%  cl.what          Text label, used in plotting routines.
%
%  sc column 1      f11.
%  sc column 2      f22.
%  sc column 3      f21 (Complex).
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%    Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Bloomfield, P. Fourier Analysis of Time Series: An Introduction.
%    2nd edition. Wiley, New York, 2000.
% 3. Nielsen J.B., Conway B.A., Halliday D.M., Perreault M-C. and Hultborn H.
%    Journal of Physiology, 569, 291-304, 2005.
%
% Type 0: [f,t,cl,sc] =  sp2_m1(0,sp1,sp2,sec_tot,samp_rate,seg_pwr,opt_str)
% Type 1: [f,t,cl,sc] =  sp2_m1(1,sp1,sp2,trig_times,duration,sec_tot,samp_rate,seg_pwr,opt_str)
% Type 2: [f,t,cl,sc] =  sp2_m1(2,sp1,sp2,trig_times,offset,seg_pts,sec_tot,samp_rate,seg_pwr,opt_str)

% Check number of output arguments
if (nargout<3)
  error(' Not enough output arguments');
end  

% Check sp_type is specified as scalar value.
if (max(size(sp_type)) ~= 1)
  error('Non scalar value for analysis type');
end

% Check for single column data
[nrow,ncol]=size(sp1);
if (ncol~=1)
  error(' Input NOT single column: sp1')
end 
[nrow,ncol]=size(sp2);
if (ncol~=1)
  error(' Input NOT single column: sp2')
end 

% Check spike data for any negative values.
if (min(sp1)<1 | min(sp2)<1)
  error(' Negative or zero spike times');
end

% Check spike data for any non integer values.
if ~isequal(round(sp1),sp1)
  error(' Non integer spike times in sp1');
end
if ~isequal(round(sp2),sp2)
  error(' Non integer spike times in sp2');
end

% Default options string.
opt_str='';

% Check and assign parameters according to analysis type.
switch (sp_type)
  case 0                                   % Type 0.
    if (nargin<6)                          % Check numbers of arguments.
      error(' Type 0 - Not enough input arguments');
    end
    sec_tot=varargin{1};     % Assign parameters.
    samp_rate=varargin{2};
    seg_pwr=varargin{3};
    if (nargin>6)
      opt_str=varargin{4};
    end
    if (max(size(sec_tot)) ~= 1)
      error(' Type 0 - Non scalar value for: sec_tot');
    end
    if (max(size(samp_rate)) ~= 1)
      error(' Type 0 - Non scalar value for: samp_rate');
    end
    if (max(size(seg_pwr)) ~= 1)
      error(' Type 0 - Non scalar value for: seg_pwr');
    end
    pts_tot=round(sec_tot*samp_rate);      % Number of samples in data set.
    seg_size=2^seg_pwr;                    % DFT segment length (S).
    seg_tot=fix(pts_tot/seg_size);         % Number of complete segments (L).
    samp_tot=seg_tot*seg_size;             % Number of samples to analyse: R=LS.
    seg_samp_start=(1:seg_size:samp_tot)'; % Start times for each segment.
    seg_samp_no(1:seg_tot,1)=seg_size;     % Fixed number of data points per segment, T = S.
    offset=0;                              % No offset values for type 0 analysis.
  case 1              % Type 1.
    if (nargin<8)     % Check numbers of arguments. 
      error(' Type 1 - Not enough input arguments');
    end  
    trig_times=varargin{1}; % Assign parameters.
    duration=varargin{2};
    sec_tot=varargin{3};
    samp_rate=varargin{4};
    seg_pwr=varargin{5};
    if (nargin>8)
      opt_str=varargin{6};
    end 
    if (max(size(sec_tot)) ~= 1)
      error(' Type 1 - Non scalar value for: sec_tot');
    end
    if (max(size(samp_rate)) ~= 1)
      error(' Type 1 - Non scalar value for: samp_rate');
    end
    if (max(size(seg_pwr)) ~= 1)
      error(' Type 1 - Non scalar value for: seg_pwr');
    end
    pts_tot=round(sec_tot*samp_rate);     % Number of samples in data set.
    seg_size=2^seg_pwr; % DFT segment length (S).
    seg_min=0.05;       % Define minimum percentage of data points allowed in each segment.
    seg_samp_min=round(seg_min*seg_size); % Convert to minimum number of data points.
    seg_tot=0;          % Counts number of segments in analysis.
    samp_tot=0;         % Counts number of samples  in analysis.
    offset=0;           % No offset values for type 1 analysis.
    [nrow,ncol]=size(trig_times);
    if (ncol~=1)
      error(' Input NOT single column: trig_times')
    end
    [nrow,ncol]=size(duration);
    if (ncol~=1)
      error(' Input NOT single column: duration')
    end
    if (min(trig_times)<1)
      error(' Type 1 - Negative or zero trig_times')
    end
    if (min(duration)<1)
      error(' Type 1 - Negative or zero duration')
    end
    if (length(trig_times) ~= length(duration))
      error(' Type 1 - Unequal numbers of trig_times, duration')
    end
    if (max(trig_times>pts_tot))
      error(' Type 1 - trig_times exceed data length')
    end
    if (max(trig_times+duration)-1>pts_tot)
      error(' Type 1 - trig_times+duration exceed data length')
    end
    for ind=1:length(trig_times)  % Loop through all trigger times.
      seg_start_offset=0;         % Offset within each block of data.
      while ((duration(ind)-seg_start_offset)>=seg_samp_min)
        seg_tot=seg_tot+1;        % Additional segment in this block.
        seg_samp_start(seg_tot,1)=trig_times(ind)+seg_start_offset;  % Start of segment.
        if ((duration(ind)-seg_start_offset)>seg_size)
          seg_samp_no(seg_tot,1)=seg_size;                        % Data for complete segment.
        else
          seg_samp_no(seg_tot,1)=duration(ind)-seg_start_offset;  % Data for part segment only.
        end 
        samp_tot=samp_tot+seg_samp_no(seg_tot,1);                 % Update number of samples.
        seg_start_offset=seg_start_offset+seg_samp_no(seg_tot,1); % Update start offset in block.
      end
    end 
  case 2                            % Type 2.
    if (nargin<9)                   % Check numbers of arguments.
      error(' Type 2 - Not enough input arguments');
    end  
    trig_times=varargin{1}; % Assign parameters.
    offset=varargin{2};
    seg_pts=varargin{3};
    sec_tot=varargin{4};
    samp_rate=varargin{5};
    seg_pwr=varargin{6};
    if (nargin>9)
      opt_str=varargin{7};
    end 
    if (max(size(sec_tot)) ~= 1)
      error(' Type 2 - Non scalar value for: sec_tot');
    end
    if (max(size(samp_rate)) ~= 1)
      error(' Type 2 - Non scalar value for: samp_rate');
    end
    if (max(size(seg_pwr)) ~= 1)
      error(' Type 2 - Non scalar value for: seg_pwr');
    end
    pts_tot=round(sec_tot*samp_rate); % Number of samples in data set.
    seg_tot=length(trig_times);       % Number of segments, L, from no of triggers.
    samp_tot=seg_tot*seg_pts;         % Number of samples to analyse.
    seg_samp_start=trig_times;        % Start times for each segment, from trigger times.
    seg_samp_no(1:seg_tot,1)=seg_pts; % Fixed number of samples per segment, T=seg_pts.
    seg_size=2^seg_pwr;               % DFT segment length (S).
    [nrow,ncol]=size(trig_times);
    if (ncol~=1)
      error(' Type 2 - Input NOT single column: trig_times')
    end
    [nrow,ncol]=size(offset);
    if (ncol~=1)
      error(' Type 2 - Input NOT single column: offset')
    end
    if (min(trig_times)<1)
      error(' Type 2 - Negative or zero trig_times')
    end
    if (max(trig_times)>pts_tot)
      error(' Type 2 - trig_times exceed data length')
    end
    if (seg_pts<1)
      error(' Type 2 - Negative or zero seg_pts')
    end
    if (seg_pts > seg_size)
      error(' Type 2 - seg_pts exceeds DFT segment size')
    end
    if (max(trig_times)+seg_pts-1>pts_tot)
      error(' Type 2 - trig_times+seg_pts exceed data length')
    end
    if (min(trig_times)+min(offset)<1)
      error(' Type 2 - Negative or zero trig_times+offset')
    end
    if (max(trig_times)+max(offset)+seg_pts-1>pts_tot)
      error(' Type 2 - trig_times+offset+seg_pts exceed data length')
    end
  otherwise
    error ([' Type error -- ',num2str(sp_type)]);  % Illegal type.
end

% Check spike times not greater than pts_tot.
if (max(sp1)>pts_tot | max(sp2)>pts_tot)
  error(' Spike times exceed specified record length')
end

% Set up 0/1 point process representations.
dat1(1:pts_tot,1)=0;
dat1(sp1)=1;
dat2(1:pts_tot,1)=0;
dat2(sp2)=1;

% Create data matrices, S rows, L columns.
rd1=zeros(seg_size,seg_tot);
rd2=zeros(seg_size,seg_tot);

% Process options.
flags.sp_type=sp_type;   % Type passed to periodogram subroutine.
if (sp_type == 2)
  flags.seg_pts=seg_pts; % Type 2: seg_pts used to define start frequency in f_out.
end 
flags.line=0;            % Set defaults - options off.
flags.inv=0;
flags.han=0;
flags.wind=0;
flags.gain=0;
flags.cov=0;
options=deblank(opt_str);
while (any(options))              % Parse individual options from string.
  [opt,options]=strtok(options);
  optarg=opt(2:length(opt));      % Determine option argument.
  switch (opt(1))
    case 'c'             % Include Periodogram COV test.
      flags.cov=1;
      if (sp_type==1)
        disp('Warning - Periodogram COV test not compatible with Type 1 analysis. Option disabled')
        flags.cov=0;
      end  
    case 'h'             % Use additional hanning smoothing on spectra.
      flags.han=1;
    case 'i'             % Channel reference invert option.
      flags.inv=1;
    case 'm'             % Mains/line frequency suppression option.
      flags.line=1;
    case 's'             % System Identification option (gain and impulse response).
      flags.gain=1;
    case 'w'             % Data window (taper) option
      i=str2num(optarg);
    if (i<1 | i>4)
      error(['error in option argument -- w',optarg]);
    end
    flags.wind=i;
    otherwise
      error (['Illegal option -- ',opt]);  % Illegal option.
  end
end

for ind_off=1:length(offset)            % Loop across all offset values.
  for ind=1:seg_tot                        % Loop across columns/segments.
    seg_pts=seg_samp_no(ind);                       % No of data points in segment.
    seg_start=seg_samp_start(ind)+offset(ind_off);  % Start sample in data vector.
    seg_stop=seg_start+seg_pts-1;                   % Stop  sample in data vector.
    dat_seg1=dat1(seg_start:seg_stop);              % Extract segment from dat1.
    dat_seg2=dat2(seg_start:seg_stop);              % Extract segment from dat2.
    md1=mean(dat_seg1);                             % Mean of segment from dat1.
    md2=mean(dat_seg2);                             % Mean of segment from dat2.
  
    rd1(1:seg_pts,ind)=dat_seg1-md1;        % Subtract mean from ch 1.
    rd2(1:seg_pts,ind)=dat_seg2-md2;        % Subtract mean from ch 2.
  end

  % Call sp2_fn2a() Weighted periodogram based spectral estimation routine.
  if (nargout>3)
    [f(:,:,ind_off),t(:,:,ind_off),cl(ind_off),sc(:,:,ind_off)] = sp2_fn2a(rd1,rd2,seg_samp_no,samp_rate,flags);
  else  
    [f(:,:,ind_off),t(:,:,ind_off),cl(ind_off)]                 = sp2_fn2a(rd1,rd2,seg_samp_no,samp_rate,flags);
  end 

  % Display No of segments & resolution
  if (sp_type ==2)
    disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),...
          ' sec,  Offset: ',num2str(offset(ind_off)),' pts,   Resolution: ',num2str(cl(ind_off).df),' Hz.']);
  else
    disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),...
          ' sec,  Resolution: ',num2str(cl(ind_off).df),' Hz.']);
  end
end

% Set additional data and type dependent elements in cl structure.
for ind_off=1:length(offset)
  switch (sp_type)
    case 0                                % Type 0.
      % Numbers of spikes analysed
      N1=length(find(sp1<cl(ind_off).samp_tot));
      N2=length(find(sp2<cl(ind_off).samp_tot));
    case 1                                % Type 1.
      % Count numbers of spikes analysed
      N1=0;
      N2=0;
      for ind_seg=1:length(trig_times)
        start=trig_times(ind_seg);
        stop=start+duration(ind_seg);
        N1=N1+length(find(sp1>=start & sp1<stop));
        N2=N2+length(find(sp2>=start & sp2<stop));
      end
    case 2                                % Type 2.
      cl(ind_off).seg_pts=seg_pts;        % Data points in segment.
      cl(ind_off).offset=offset(ind_off); % Offset.
      % Count numbers of spikes analysed
      N1=0;
      N2=0;
      for ind_seg=1:length(trig_times)
        start=trig_times(ind_seg)+offset(ind_off);
        stop=start+seg_pts;
        N1=N1+length(find(sp1>=start & sp1<stop));
        N2=N2+length(find(sp2>=start & sp2<stop));
      end
  end
  if (flags.inv)      % Input and output channels inverted with 'i' option.
    cl(ind_off).N2=N1;                             % N2, No of events analysed from sp1.
    cl(ind_off).N1=N2;                             % N1, No of events analysed from sp2.
    cl(ind_off).P2=N1/(cl(ind_off).samp_tot*2*pi); % P2, mean intensity of events from sp1.
    cl(ind_off).P1=N2/(cl(ind_off).samp_tot*2*pi); % P1, mean intensity of events from sp2.
  else
    cl(ind_off).N1=N1;                             % N1, No of events analysed from sp1.
    cl(ind_off).N2=N2;                             % N2, No of events analysed from sp2.
    cl(ind_off).P1=N1/(cl(ind_off).samp_tot*2*pi); % P1, mean intensity of events from sp1.
    cl(ind_off).P2=N2/(cl(ind_off).samp_tot*2*pi); % P2, mean intensity of events from sp2.
  end
  cl(ind_off).opt_str=opt_str;                     % Copy of options string.
  cl(ind_off).what='';                             % Field for plot label.
end

function u = simNoiseIntr(R, m)
% simulate the instrinsic noise - based on BAA_plotModelFit

% calls the background noise initializer and normalize by proper constant
u = innovate_timeseries(R,m);
u{1} = u{1}.*sqrt(R.IntP.dt);
end
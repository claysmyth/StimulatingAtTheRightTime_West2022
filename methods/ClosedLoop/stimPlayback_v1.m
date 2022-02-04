function [uexs,R,phi] = stimPlayback_v1(~,R,~,~,~,~,phi)
load([R.rootn 'data\phaseStimSave\stim_tmp_' sprintf('%3.f',1000*R.IntP.phaseStim.phaseshift)],'uexs');
uexs = uexs';
uexs = flipud(uexs);
R.IntP.phaseStim.switch = 0; % turn off the stim update
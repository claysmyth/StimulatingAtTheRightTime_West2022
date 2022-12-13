    % now test equality for debugging purposes
    try
        for i = 1:numEpoch
            assert(allclose(xsims_stim_Cond_SPrev{i}, vec_xsims_SPrevComp{i}));
        end
    catch
        t1 = 1;
    end

    t1 = xsims_gl_stim_Cond{i};
    t2 = vec_xsims_glComp{i};
    t3 = uexsS_Cond{i};

    t1 = resample(t1', 1, 10);
    t2 = resample(t2', 1, 10);

    hpFilt = designfilt('bandpassfir','FilterOrder', 20, ...
        'CutoffFrequency1', 13,'CutoffFrequency2', 20, ...
         'SampleRate',200);
    hpFilt2 = designfilt('highpassfir','FilterOrder', 20, ...
         'CutoffFrequency', 5,'PassbandRipple', 1, ...
         'SampleRate',200);

    figure; title('Original LFP');
    for i = 1:6
        subplot(7, 1, i); hold on;
        plot(t1(:, i), 'DisplayName', 'withStim'); 
        plot(t2(:, i), 'DisplayName', 'Original');
        legend('boxoff');
    end
    subplot(717); plot(t3(:, 4));

    t1Beta = filtfilt(hpFilt, t1)';
    t2Beta = filtfilt(hpFilt, t2)';

    t1Theta = filtfilt(hpFilt2, t1)';
    t2Theta = filtfilt(hpFilt2, t2)';

    figure; 
    for i = 1:6
        subplot(7, 1, i); hold on;
        plot(t1Beta(i, :), 'DisplayName', 'withStim'); 
        plot(t2Beta(i, :), 'DisplayName', 'Original'); ylim([-10, 10])
        legend('boxoff');
    end
    subplot(717); plot(t3(:, R.IntP.phaseStim.sensStm(2))); title('Everything above theta');

    figure; 
    for i = 1:6
        subplot(7, 1, i); hold on;
        plot(t1Theta(i, :), 'DisplayName', 'withStim'); 
        plot(t2Theta(i, :), 'DisplayName', 'Original'); ylim([-10, 10])
        legend('boxoff');
    end
    subplot(717); plot(t3(:, R.IntP.phaseStim.sensStm(2))); title('Everything above beta');

    % perform hiltert
    t1BetaAmp = abs(hilbert(t1Beta'))';
    t2BetaAmp = abs(hilbert(t2Beta'))';
    figure; 
    for i = 1:6
        subplot(7, 1, i); hold on;
        plot(t1BetaAmp(i, :), 'DisplayName', 'withStim'); 
        plot(t2BetaAmp(i, :), 'DisplayName', 'Original'); ylim([-10, 10])
        legend('boxoff');
    end
    subplot(717); plot(t3(:, R.IntP.phaseStim.sensStm(2))); title('Everything above beta hilbert');
    t1 = 1;
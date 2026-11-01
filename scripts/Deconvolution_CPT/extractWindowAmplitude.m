function amp = extractWindowAmplitude(EEG, window_ms, channels)
    chanIdx = find(ismember({EEG.chanlocs.labels}, channels));
    if isempty(chanIdx)
        error('None of the requested ROI channels were found in this subject''s data.');
    end
    timeIdx = EEG.times >= window_ms(1) & EEG.times <= window_ms(2);
    amp = nan(EEG.trials, 1);
    for e = 1:EEG.trials
        amp(e) = mean(mean(EEG.data(chanIdx, timeIdx, e), 1), 2);
    end
end
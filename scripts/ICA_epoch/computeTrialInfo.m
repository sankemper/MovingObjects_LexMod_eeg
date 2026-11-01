function T = computeTrialInfo(EEGepoched, posLabel, subjectID, windows, roi)
    n = EEGepoched.trials;
    epochIdx = (1:n)';
    item = nan(n,1); mention = nan(n,1); modified = false(n,1);
    pattern = cell(n,1); informative = false(n,1); block = nan(n,1);
    for e = 1:n
        lat = EEGepoched.epoch(e).eventlatency;
        if iscell(lat), lat = cell2mat(lat); end
        [~, targetIdx] = min(abs(lat));
        item(e)        = EEGepoched.epoch(e).eventitem{targetIdx};
        mention(e)     = EEGepoched.epoch(e).eventmention{targetIdx};
        modified(e)    = EEGepoched.epoch(e).eventmodified{targetIdx};
        pattern{e}     = EEGepoched.epoch(e).eventpattern{targetIdx};
        informative(e) = EEGepoched.epoch(e).eventinformative{targetIdx};
        block(e)       = EEGepoched.epoch(e).eventblock{targetIdx};
    end

    T = table(repmat({subjectID}, n, 1), repmat({posLabel}, n, 1), epochIdx, ...
        item, mention, modified, pattern, informative, block, ...
        'VariableNames', {'subject','position','epoch','item','mention','modified', ...
        'pattern','informative','block'});

    % one column per window; each window name needs a matching ROI name
    for nm = fieldnames(windows)'
        T.(nm{1}) = extractWindowAmplitude(EEGepoched, windows.(nm{1}), roi.(nm{1}));
    end
end
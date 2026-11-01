function chanlocsOut = merge_locations_by_label(chanlocsIn, chanlocsSource)
    locFields = {'theta','radius','X','Y','Z','sph_theta','sph_phi','sph_radius'};
    presentFields = locFields(isfield(chanlocsSource, locFields));
    if isempty(presentFields)
        error('merge_locations_by_label:noGeomFields', ...
            ['chanlocsSource has none of the expected geometry fields. ' ...
             'Run fieldnames(chanlocsSource) and check what loadbvef actually returned.']);
    end

    srcLabels = lower({chanlocsSource.labels});
    chanlocsOut = chanlocsIn;
    nMatched = 0;

    for i = 1:numel(chanlocsOut)
        idx = find(strcmpi(chanlocsOut(i).labels, srcLabels), 1);
        if ~isempty(idx)
            for f = 1:numel(presentFields)
                fn = presentFields{f};
                chanlocsOut(i).(fn) = chanlocsSource(idx).(fn);
            end
            nMatched = nMatched + 1;
        end
    end

    fprintf('merge_locations_by_label: matched %d / %d channels to a scalp position\n', ...
        nMatched, numel(chanlocsOut));
    unmatched = {chanlocsOut(~ismember(lower({chanlocsOut.labels}), srcLabels)).labels};
    if ~isempty(unmatched)
        fprintf('  no position found for: %s\n', strjoin(unmatched, ', '));
    end
end
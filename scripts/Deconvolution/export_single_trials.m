function T = export_single_trials(EEG, subjectID, badSamples, epochWin, chanIdx, binMs)
fs   = EEG.srate;
uf   = EEG.unfold;
nPred = size(uf.X, 2);
nLag  = size(uf.beta_dc, 2);
nCh   = numel(chanIdx);

assert(size(uf.Xdc,2) == nPred*nLag, 'Xdc columns (%d) ~= predictors x lags (%d).', size(uf.Xdc,2), nPred*nLag);
assert(size(uf.Xdc,1) == EEG.pnts,   'Xdc rows (%d) ~= EEG.pnts (%d).', size(uf.Xdc,1), EEG.pnts);
assert(size(uf.X,1)   == numel(EEG.event), 'X rows (%d) ~= number of EEG.event (%d).', size(uf.X,1), numel(EEG.event));

% event type name per predictor
names = cell(1, numel(uf.eventtypes));
for k = 1:numel(uf.eventtypes)
    nm = uf.eventtypes{k};
    while iscell(nm), nm = nm{1}; end
    names{k} = char(nm);
end
predEv = names(uf.cols2eventtypes(:)');            % 1 x nPred
colEv  = repelem(predEv, nLag);                    % 1 x nPred*nLag (predictor-major, lags contiguous)

% betas as channels x (nPred*nLag), same ordering as Xdc columns
B = reshape(uf.beta_dc(chanIdx, :, :), nCh, nLag*nPred);

evTypes = {'DET1_crit','ADJ1_crit','N1_crit'};
labels  = string({EEG.chanlocs(chanIdx).labels})';
evAll   = EEG.event;
evTyp   = string({evAll.type});
s0 = round(epochWin(1)*fs);  s1 = round(epochWin(2)*fs);
tAll = (s0:s1) / fs;
bl   = tAll >= -0.1 & tAll < 0;
parts = {};

for e = 1:numel(evTypes)
    other = ~strcmp(colEv, evTypes{e});
    % modeled signal from all OTHER event types (samples x channels -> channels x samples)
    predOther = (uf.Xdc(:, other) * B(:, other)')';
    cleaned   = EEG.data(chanIdx, :) - predOther;
    clear predOther

    ev  = evAll(evTyp == evTypes{e});
    lat = round([ev.latency]);
    keep = lat + s0 >= 1 & lat + s1 <= EEG.pnts;

    X = nan(nCh, numel(tAll), numel(lat));
    for k = find(keep)
        idx = lat(k) + (s0:s1);
        if any(badSamples(idx)), keep(k) = false; continue; end
        ep = cleaned(:, idx);
        X(:,:,k) = ep - mean(ep(:, bl), 2);               % per-trial baseline
    end
    X = X(:,:,keep);  evK = ev(keep);  nTr = numel(evK);
    fprintf('%s %s: %d of %d trials kept\n', subjectID, evTypes{e}, nTr, numel(ev));

    [Xb, tb] = bin_time(X, tAll, binMs);
    nTb = numel(tb);
    [chI, tI, cI] = ndgrid(1:nCh, 1:nTb, 1:nTr);

    T = table( ...
    repmat(string(subjectID), numel(chI), 1), ...
    repmat(string(evTypes{e}), numel(chI), 1), ...
    cI(:), evfield(evK,'criticalTrialID',cI), evfield(evK,'pattern',cI), ...
    evfield(evK,'mention',cI), evfield(evK,'modified',cI), evfield(evK,'informative',cI), ...
    labels(chI(:)), tb(tI(:))' * 1000, Xb(:), ...
    'VariableNames', {'subject','event','trial_idx','criticalTrialID','pattern', ...
    'mention','modified','informative','channel','time','amplitude'});
    parts{end+1} = T; %#ok<AGROW>
end
T = vertcat(parts{:});
end

function v = evfield(ev, f, cI)
if ~isfield(ev, f)
    warning('EEG.event has no field "%s"; filling with NA.', f);
    v = repmat("NA", numel(cI), 1); return
end
c = {ev.(f)}';
if all(cellfun(@(x) isnumeric(x) && isscalar(x), c))
    vals = cell2mat(c);
else
    c(cellfun(@isempty, c)) = {"NA"};
    vals = cellfun(@(x) string(x), c);
end
v = vals(cI(:));
end

function [out, tb] = bin_time(data, times, binMs)
if binMs <= 0, out = data; tb = times; return, end
bi = floor((times - times(1)) * 1000 / binMs) + 1;
nb = max(bi);
out = zeros(size(data,1), nb, size(data,3));
tb  = zeros(1, nb);
for b = 1:nb
    m = bi == b;
    out(:, b, :) = mean(data(:, m, :), 2);
    tb(b) = mean(times(m));
end
end
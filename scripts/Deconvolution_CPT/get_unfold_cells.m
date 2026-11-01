function cells = get_unfold_cells(ufresult, eventName, predNames, levels, baselineWin)
% GET_UNFOLD_CELLS  Modeled waveform for every design cell of one event type,
% computed directly from the treatment-coded betas (no uf_predict* needed).
%
%   cells = get_unfold_cells(ufresult, 'N1_crit', {'mention','modified','informative'}, ...
%                            {[1 2],[0 1],[0 1]}, [-0.1 0])
%
% A cell's waveform = intercept + every beta whose dummy terms are all "on"
% for that cell (so interaction betas are included automatically).
%
% Output: cells.data (chan x time x cell), cells.design (table), cells.times, cells.event

if nargin < 5, baselineWin = []; end

% ---- parameters belonging to this event --------------------------------
nPar = numel(ufresult.param);
pIdx = [];
for i = 1:nPar
    ev = ufresult.param(i).event;
    if iscell(ev), ev = ev{1}; end
    if strcmp(char(ev), eventName), pIdx(end+1) = i; end %#ok<AGROW>
end
assert(~isempty(pIdx), 'No parameters found for event %s.', eventName);
assert(size(ufresult.beta, 3) == nPar, 'beta 3rd dim (%d) ~= numel(param) (%d).', ...
       size(ufresult.beta, 3), nPar);

names = {ufresult.param(pIdx).name};
fprintf('%s coefficients: %s\n', eventName, strjoin(names, ' | '));

% ---- parse each coefficient into {variable, level} dummy terms ---------
terms = cell(1, numel(pIdx));      % each: struct array of var/level (empty = intercept)
for k = 1:numel(pIdx)
    nm = names{k};
    if strcmpi(ufresult.param(pIdx(k)).type, 'intercept') || contains(lower(nm), 'intercept')
        terms{k} = struct('var', {}, 'level', {});
        continue
    end
    parts = strsplit(nm, ':');
    t = struct('var', cell(1, numel(parts)), 'level', cell(1, numel(parts)));
    for q = 1:numel(parts)
        pq  = regexprep(strtrim(parts{q}), '^\d+_', '');     % event-index prefix, e.g. '2_'
        pq  = regexprep(pq, 'cat\(([^)]*)\)', '$1');
        tok = regexp(pq, '^(.+?)[_=]\(?([^_=()]+)\)?$', 'tokens', 'once');
        if isempty(tok) || ~ismember(tok{1}, predNames)
            error(['get_unfold_cells: cannot parse coefficient name "%s". ' ...
                   'Paste the coefficient list printed above.'], nm);
        end
        t(q).var   = tok{1};
        t(q).level = level2num(tok{2});
    end
    terms{k} = t;
end

% ---- all cells (full factorial, first predictor varies fastest) --------
grids = cell(1, numel(levels));
[grids{:}] = ndgrid(levels{:});
nCell = numel(grids{1});
design = table();
for k = 1:numel(predNames)
    design.(predNames{k}) = double(grids{k}(:));
end

% ---- sum betas per cell -------------------------------------------------
B = ufresult.beta(:, :, pIdx);                  % chan x time x par
data = zeros(size(B,1), size(B,2), nCell);
for c = 1:nCell
    acc = zeros(size(B,1), size(B,2));
    for k = 1:numel(pIdx)
        t = terms{k};
        on = true;
        for q = 1:numel(t)
            on = on && design.(t(q).var)(c) == t(q).level;
        end
        if on, acc = acc + B(:, :, k); end
    end
    data(:, :, c) = acc;
end

times = ufresult.times(:)';
if ~isempty(baselineWin)
    b = times >= baselineWin(1) & times <= baselineWin(2);
    data = data - mean(data(:, b, :), 2);
end

cells.data   = data;
cells.design = design;
cells.times  = times;
cells.event  = eventName;
end

function x = level2num(s)
s = lower(strtrim(s));
if strcmp(s, 'true'),  x = 1;
elseif strcmp(s, 'false'), x = 0;
else, x = str2double(s);
end
end
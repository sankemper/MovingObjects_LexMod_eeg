function cells = save_subject_cells(ufresult, subjectID, outFolder, baselineWin)
% Baseline once, from the DET1 pre-onset window (clean fixation period), applied
% to all three events. Only intercepts are shifted: every cell waveform is
% intercept + effect terms, so this shifts each cell by the same constant
% and leaves the condition contrasts untouched.
if nargin < 4, baselineWin = []; end
if ~isempty(baselineWin)
    ufresult = baseline_from_det1(ufresult, baselineWin);
end

cells.DET1_crit = get_unfold_cells(ufresult, 'DET1_crit', ...
    {'mention','modified','informative'}, {[1 2],[0 1],[0 1]}, []);
% adj1: 'modified' was dropped from the formula (constant)
cells.ADJ1_crit = get_unfold_cells(ufresult, 'ADJ1_crit', ...
    {'mention','informative'}, {[1 2],[0 1]}, []);
cells.N1_crit = get_unfold_cells(ufresult, 'N1_crit', ...
    {'mention','modified','informative'}, {[1 2],[0 1],[0 1]}, []);
save(fullfile(outFolder, sprintf('%s_unfold_cells.mat', subjectID)), 'cells');
end

function ufresult = baseline_from_det1(ufresult, win)
% Subtract one baseline per channel, taken from the DET1 intercept's pre-onset
% window, from ALL intercepts. Each cell waveform is intercept + effect terms,
% so every cell of every event shifts by the same constant and the condition
% contrasts stay untouched.
nPar = numel(ufresult.param);
isInt = false(1, nPar);  isDet = false(1, nPar);
for i = 1:nPar
    ev = ufresult.param(i).event; if iscell(ev), ev = ev{1}; end
    nm = ufresult.param(i).name;  if iscell(nm), nm = nm{1}; end
    ty = ufresult.param(i).type;  if iscell(ty), ty = ty{1}; end
    isInt(i) = strcmpi(char(ty), 'intercept') || contains(lower(char(nm)), 'intercept');
    isDet(i) = strcmp(char(ev), 'DET1_crit');
end
detInt = find(isInt & isDet);
assert(numel(detInt) == 1, 'Expected exactly one DET1_crit intercept, found %d.', numel(detInt));

tIdx = ufresult.times >= win(1) & ufresult.times <= win(2);
base = mean(ufresult.beta(:, tIdx, detInt), 2);          % channels x 1
for p = find(isInt)
    ufresult.beta(:, :, p) = ufresult.beta(:, :, p) - base;
end
end
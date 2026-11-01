function T = export_unfold_gamm(outFolder, subjectNums, binMs, outFile)
if nargin < 3, binMs = 20; end
if nargin < 4, outFile = 'unfold_gamm_long.csv'; end
tmp = pop_loadset('filename', sprintf('P%d_ica_processed.set', subjectNums(1)), ...
                  'filepath', outFolder, 'loadmode', 'info');
chanLabels = {tmp.chanlocs.labels}';

events = {'DET1_crit','ADJ1_crit','N1_crit'};
parts  = {};

for s = subjectNums
    subjectID = sprintf('P%d', s);
    S = load(fullfile(outFolder, sprintf('%s_unfold_cells.mat', subjectID)));

    for e = 1:numel(events)
        c = S.cells.(events{e});
        [nCh, ~, nC] = size(c.data);
        assert(nCh == numel(chanLabels), 'Channel count mismatch (%d vs %d).', nCh, numel(chanLabels));

        [data, tb] = bin_time(c.data, c.times, binMs);
        nT = numel(tb);

        design = c.design;
        if ~ismember('modified', design.Properties.VariableNames)
            design.modified = ones(height(design), 1);   % adj1: always modified
        end

        [chI, tI, cI] = ndgrid(1:nCh, 1:nT, 1:nC);
        n = numel(chI);
        T = table( ...
            repmat(string(subjectID), n, 1), ...
            repmat(string(events{e}),  n, 1), ...
            design.mention(cI(:)), ...
            design.modified(cI(:)), ...
            design.informative(cI(:)), ...
            string(chanLabels(chI(:))), ...
            tb(tI(:))' * 1000, ...                       % ms
            data(:), ...
            'VariableNames', {'subject','event','mention','modified', ...
                              'informative','channel','time','amplitude'});
        parts{end+1} = T; %#ok<AGROW>
    end
    fprintf('%s exported\n', subjectID);
end

T = vertcat(parts{:});
writetable(T, fullfile(outFolder, outFile));
fprintf('Wrote %d rows to %s\n', height(T), fullfile(outFolder, outFile));
end

function [out, tb] = bin_time(data, times, binMs)
if binMs <= 0
    out = data; tb = times; return
end
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
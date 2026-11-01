function [clusterTbl, testTbl] = summarize_unfold_clusters(outFolder, elec)
% SUMMARIZE_UNFOLD_CLUSTERS  Log and plot every unfold_stat_*.mat cluster test.
%
%   [clusterTbl, testTbl] = summarize_unfold_clusters(outFolder, elec);
%
% elec = FieldTrip electrode structure (label, elecpos, chanpos, unit),
% the same one you pass to unfold_cluster_test.
%
% Outputs (all in <outFolder>/ClusterSummary):
%   cluster_log_all_clusters.csv  one row per cluster (every cluster found, significant or not)
%   cluster_log_tests.csv         one row per test (smallest p, FDR within family)
%   cluster_summary.txt           readable log of everything
%   plots/                        topography series + channel-by-time map per test
%
% P-VALUES: unfold_cluster_test uses tail = 0 and alpha = 0.025, so the
% p-values FieldTrip returns are one-sided per direction. This script
% reports p_two_sided = min(1, 2 * p) (compare to 0.05). If you later add
% cfg.correcttail = 'prob' inside unfold_cluster_test, this is detected
% automatically and the doubling is skipped.
%
% FAMILIES: tests whose name contains '_x_' (interactions, including the
% three-way) form one family; all other tests (main effects) form another.
% FDR (Benjamini-Hochberg) is applied within each family. Tests with no
% cluster at all count as p = 1 in their family, so the family size equals
% the number of tests you actually ran.

outDir  = fullfile(outFolder, 'ClusterSummary');
plotDir = fullfile(outDir, 'plots');
if ~exist(plotDir, 'dir'), mkdir(plotDir); end

files = dir(fullfile(outFolder, 'unfold_stat_*.mat'));
if isempty(files)
    error('No unfold_stat_*.mat files found in %s', outFolder);
end

clusterRows = {};
testRows    = {};
layout      = [];

for i = 1:numel(files)
    tok = regexp(files(i).name, '^unfold_stat_(DET1_crit|ADJ1_crit|N1_crit)_(.+)\.mat$', 'tokens', 'once');
    if isempty(tok)
        warning('Skipping %s (name does not match unfold_stat_<event>_<label>.mat)', files(i).name);
        continue;
    end
    eventName = tok{1};
    label     = tok{2};

    S    = load(fullfile(outFolder, files(i).name));
    stat = S.stat;

    % layout built once, from the channels in the first stat file
    if isempty(layout)
        [tf, loc] = ismember(stat.label, elec.label);
        assert(all(tf), 'Not in elec: %s', strjoin(stat.label(~tf)', ', '));
        elecSub = struct('label', {elec.label(loc)}, 'elecpos', elec.elecpos(loc,:), ...
                         'chanpos', elec.chanpos(loc,:), 'unit', elec.unit);
        layout = ft_prepare_layout(struct('elec', elecSub));
    end

    % number of subjects, from the stored design (2 rows x 2n columns)
    nSubj = NaN;
    if isfield(stat, 'cfg') && isfield(stat.cfg, 'design') && ~isempty(stat.cfg.design)
        nSubj = size(stat.cfg.design, 2) / 2;
    end

    % one-sided-per-tail -> two-sided, unless FieldTrip already did it
    alreadyTwoSided = isfield(stat, 'cfg') && isfield(stat.cfg, 'correcttail') ...
        && ischar(stat.cfg.correcttail) && strcmp(stat.cfg.correcttail, 'prob');
    factor = 2;
    if alreadyTwoSided, factor = 1; end

    pAll = [];
    nPos = 0;
    nNeg = 0;
    for pol = {'pos', 'neg'}
        fName = [pol{1} 'clusters'];
        lName = [pol{1} 'clusterslabelmat'];
        if ~isfield(stat, fName) || isempty(stat.(fName)) || ~isfield(stat, lName)
            continue;
        end
        for k = 1:numel(stat.(fName))
            m    = stat.(lName) == k;
            pRaw = stat.(fName)(k).prob;
            pTwo = min(1, pRaw * factor);
            pAll(end+1) = pTwo; %#ok<AGROW>
            if strcmp(pol{1}, 'pos'), nPos = nPos + 1; else, nNeg = nNeg + 1; end

            if any(m(:))
                tv = stat.stat;
                tv(~m) = NaN;
                [~, ix] = max(abs(tv(:)));
                [ci, ti] = ind2sub(size(tv), ix);
                timesIn = stat.time(any(m, 1));
                chansIn = reshape(stat.label(any(m, 2)), 1, []);
                peakT = stat.stat(ci, ti);
                peakCh = stat.label{ci};
                peakMs = stat.time(ti) * 1000;
                startMs = min(timesIn) * 1000;
                endMs = max(timesIn) * 1000;
                nCh = sum(any(m, 2));
                chStr = strjoin(chansIn, ' ');
            else
                peakT = NaN; peakCh = ''; peakMs = NaN;
                startMs = NaN; endMs = NaN; nCh = 0; chStr = '';
            end

            clusterRows{end+1} = struct('test', label, 'event', eventName, ...
                'polarity', pol{1}, 'cluster', k, ...
                'p_raw', pRaw, 'p_two_sided', pTwo, ...
                'cluster_stat', stat.(fName)(k).clusterstat, ...
                'peak_t', peakT, 'peak_channel', peakCh, 'peak_time_ms', peakMs, ...
                'start_ms', startMs, 'end_ms', endMs, ...
                'n_channels', nCh, 'channels', chStr); %#ok<AGROW>
        end
    end

    if isempty(pAll)
        pMin = NaN;  pUsed = 1;      % no cluster found -> counts as p = 1 in the family
    else
        pMin = min(pAll);  pUsed = pMin;
    end
    family = 'main';
    if contains(label, '_x_'), family = 'interaction'; end

    testRows{end+1} = struct('test', label, 'event', eventName, 'family', family, ...
        'n_subjects', nSubj, 'n_pos_clusters', nPos, 'n_neg_clusters', nNeg, ...
        'min_p_two_sided', pMin, 'p_used', pUsed, 'file', files(i).name); %#ok<AGROW>

    % ---- plots ----
    try
        plot_topo_series(stat, S, layout, label, eventName, nSubj, plotDir);
        plot_raster(stat, elec, label, eventName, nSubj, plotDir);
    catch ME
        warning('Plotting failed for %s | %s: %s', label, eventName, ME.message);
    end
    fprintf('processed %s | %s\n', label, eventName);
end

% ---- tables ----
testTbl = struct2table(vertcat(testRows{:}));
testTbl.p_FDR = nan(height(testTbl), 1);
for fam = {'main', 'interaction'}
    idx = strcmp(testTbl.family, fam{1});
    if any(idx)
        testTbl.p_FDR(idx) = fdr_bh(testTbl.p_used(idx));
    end
end
testTbl.sig_uncorrected = testTbl.min_p_two_sided < 0.05;
testTbl.sig_FDR         = testTbl.p_FDR < 0.05;
testTbl = sortrows(testTbl, {'event', 'family', 'test'});

hasClusters = ~isempty(clusterRows);
if hasClusters
    clusterTbl = struct2table(vertcat(clusterRows{:}));
    clusterTbl.sig_uncorrected = clusterTbl.p_two_sided < 0.05;
    clusterTbl = sortrows(clusterTbl, {'event', 'test', 'p_two_sided'});
else
    clusterTbl = table();
end

writetable(testTbl, fullfile(outDir, 'cluster_log_tests.csv'));
if hasClusters
    writetable(clusterTbl, fullfile(outDir, 'cluster_log_all_clusters.csv'));
end

% ---- readable text log (console + file) ----
fid = fopen(fullfile(outDir, 'cluster_summary.txt'), 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
tee(fid, 'Cluster-based permutation summary (%s)\n', datestr(now));
tee(fid, 'p_two_sided = one-sided-per-tail p * 2; FDR within family (main / interaction)\n');
for r = 1:height(testTbl)
    tee(fid, '\n=== %s | %s | family: %s | n = %g ===\n', testTbl.test{r}, testTbl.event{r}, ...
        testTbl.family{r}, testTbl.n_subjects(r));
    if hasClusters
        sel = strcmp(clusterTbl.test, testTbl.test{r}) & strcmp(clusterTbl.event, testTbl.event{r});
        rows = find(sel);
    else
        rows = [];
    end
    if isempty(rows)
        tee(fid, '  no clusters found\n');
    end
    for q = rows'
        tee(fid, '  %s cluster %d: p = %.4f (two-sided), t_peak = %.2f at %s, %.0f-%.0f ms, %d channels: %s\n', ...
            clusterTbl.polarity{q}, clusterTbl.cluster(q), clusterTbl.p_two_sided(q), ...
            clusterTbl.peak_t(q), clusterTbl.peak_channel{q}, clusterTbl.start_ms(q), ...
            clusterTbl.end_ms(q), clusterTbl.n_channels(q), clusterTbl.channels{q});
    end
    tee(fid, '  smallest p = %.4f | FDR-adjusted (within %s family) = %.4f\n', ...
        testTbl.min_p_two_sided(r), testTbl.family{r}, testTbl.p_FDR(r));
end

disp(testTbl(:, {'test', 'event', 'family', 'n_subjects', 'min_p_two_sided', 'p_FDR', 'sig_FDR'}));
fprintf('\nLogs and plots saved to: %s\n', outDir);
end

%% ========================================================================
%% Plot 1: topographies over time in 50 ms steps, significant channels marked
%% ========================================================================
function plot_topo_series(stat, S, layout, label, eventName, nSubj, plotDir)
    tmax  = max(stat.time);
    edges = 0:0.05:tmax;
    if numel(edges) < 2, return; end
    nWin  = numel(edges) - 1;

    % plot the amplitude difference wave if the stat file contains it,
    % otherwise the t-values from the statistics
    if isfield(S, 'GA_diff')
        data = S.GA_diff;  param = 'avg';   unitLabel = '\muV';
        tIdx = data.time >= 0 & data.time <= tmax;
        vals = data.avg(:, tIdx);
    else
        data = stat;       param = 'stat';  unitLabel = 't';
        vals = stat.stat;
    end
    vals = vals(~isnan(vals));
    zl = max(abs(vals));
    if isempty(zl) || zl == 0, zl = 1; end

    nCols = 5;
    nRows = ceil(nWin / nCols);
    fig = figure('Visible', 'off', 'Position', [50 50 1500 280 * nRows]);
    for k = 1:nWin
        subplot(nRows, nCols, k);
        cfgT = [];
        cfgT.figure      = 'gca';
        cfgT.parameter   = param;
        cfgT.xlim        = [edges(k) edges(k+1)];
        cfgT.zlim        = [-zl zl];
        cfgT.layout      = layout;
        cfgT.colormap    = diverging_cmap();
        cfgT.comment     = 'no';
        cfgT.marker      = 'off';
        cfgT.interactive = 'no';
        cfgT.colorbar    = 'no';

        inWin = stat.time >= edges(k) & stat.time <= edges(k+1);
        sig = any(stat.mask(:, inWin), 2);
        if any(sig)
            cfgT.highlight        = 'on';
            cfgT.highlightchannel = stat.label(sig);
            cfgT.highlightsymbol  = '.';
            cfgT.highlightsize    = 14;
        else
            cfgT.highlight = 'off';
        end
        ft_topoplotER(cfgT, data);
        title(sprintf('%d-%d ms', round(edges(k)*1000), round(edges(k+1)*1000)), 'FontSize', 9);
    end
    cb = colorbar('Position', [0.94 0.3 0.01 0.4]);
    cb.Label.String = unitLabel;
    try
        sgtitle(sprintf('%s | %s | n = %g | dots = channels in a significant cluster', ...
            label, eventName, nSubj), 'Interpreter', 'none');
    catch
    end
    save_and_close(fig, plotDir, sprintf('%s_%s_topo', eventName, label));
end

%% ========================================================================
%% Plot 2: channel x time t-map, significant clusters outlined
%% ========================================================================
function plot_raster(stat, elec, label, eventName, nSubj, plotDir)
    t = stat.time * 1000;

    % order channels front -> back using the electrode x coordinate
    [~, loc] = ismember(stat.label, elec.label);
    [~, order] = sort(elec.elecpos(loc, 1), 'descend');

    tv = stat.stat(order, :);
    mk = stat.mask(order, :);
    nCh = numel(order);

    fig = figure('Visible', 'off', 'Position', [100 100 1000 560]);
    imagesc(t, 1:nCh, tv);
    colormap(gca, diverging_cmap());
    zl = max(abs(tv(~isnan(tv))));
    if isempty(zl) || zl == 0, zl = 1; end
    set(gca, 'CLim', [-zl zl]);
    hold on;
    if any(mk(:))
        contour(t, 1:nCh, double(mk), [0.5 0.5], 'k', 'LineWidth', 1.5);
    end
    set(gca, 'YTick', 1:nCh, 'YTickLabel', stat.label(order));
    xlabel('Time (ms)');
    ylabel('Channel (front to back)');
    cb = colorbar;
    cb.Label.String = 't-value';
    title(sprintf('%s | %s | n = %g | black outline = significant cluster', ...
        label, eventName, nSubj), 'Interpreter', 'none');
    save_and_close(fig, plotDir, sprintf('%s_%s_raster', eventName, label));
end

%% ========================================================================
%% Helpers
%% ========================================================================
function save_and_close(fig, plotDir, fileName)
    pngPath = fullfile(plotDir, [fileName '.png']);
    try
        exportgraphics(fig, pngPath, 'Resolution', 300);   % R2020a or newer
    catch
        saveas(fig, pngPath);
    end
    savefig(fig, fullfile(plotDir, [fileName '.fig']));
    close(fig);
end

function cmap = diverging_cmap()
    n = 128;
    blue2white = [linspace(0.1, 1, n)', linspace(0.3, 1, n)', linspace(0.8, 1, n)'];
    white2red  = [linspace(1, 0.8, n)', linspace(1, 0.1, n)', linspace(1, 0.1, n)'];
    cmap = [blue2white; white2red];
end

function q = fdr_bh(p)
    % Benjamini-Hochberg adjusted p-values
    p = p(:);
    n = numel(p);
    [ps, si] = sort(p);
    qs = ps .* n ./ (1:n)';
    qs = flipud(cummin(flipud(qs)));
    qs = min(qs, 1);
    q = zeros(size(p));
    q(si) = qs;
end

function tee(fid, fmt, varargin)
    fprintf(1, fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end
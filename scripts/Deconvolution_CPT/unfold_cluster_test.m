function stat = unfold_cluster_test(outFolder, subjectNums, eventName, ...
                                    factA, lvlA, label, elec, factB, lvlB, factC, lvlC)
% UNFOLD_CLUSTER_TEST  Dependent-samples cluster permutation test on modeled waveforms.
%
% Main effect (A level 1 vs level 2, averaged over all other factors):
%   unfold_cluster_test(f, subj, 'N1_crit', 'modified', [1 0], 'Modified_Noun', elec)
%
% Interaction (A x B, difference of differences):
%   unfold_cluster_test(f, subj, 'N1_crit', 'modified', [1 0], 'Modified_x_Mention', ...
%                       elec, 'mention', [1 2])
%
% Three-way interaction (A x B x C): difference of the A x B interaction
% between the two levels of C:
%   unfold_cluster_test(f, subj, 'N1_crit', 'modified', [1 0], 'Modified_x_Mention_x_Informative', ...
%                       elec, 'mention', [1 2], 'informative', [1 0])
%
% Match alpha / numrandomization / minnbchan below to your existing
% run_*_cluster_test functions so results are comparable.

isInteraction = nargin >= 9;
isThreeWay    = nargin >= 11;
n = numel(subjectNums);
dataA = cell(1, n);
dataB = cell(1, n);

tmp = pop_loadset('filename', sprintf('P%d_ica_processed.set', subjectNums(1)), ...
                  'filepath', outFolder, 'loadmode', 'info');
chanLabels = {tmp.chanlocs.labels};
[tf, loc] = ismember(chanLabels, elec.label);
fprintf('[1] channel labels loaded\n');
assert(all(tf), 'Not in elec: %s', strjoin(chanLabels(~tf), ', '));
elecData = struct('label', {elec.label(loc)}, 'elecpos', elec.elecpos(loc,:), ...
                  'chanpos', elec.chanpos(loc,:), 'unit', elec.unit);

for i = 1:n
    S = load(fullfile(outFolder, sprintf('P%d_unfold_cells.mat', subjectNums(i))));
    fprintf('[2] loaded subject %d/%d\n', i, n);
    c = S.cells.(eventName);
    d = c.design;
    nCh = size(c.data, 1);
    assert(nCh == numel(elecData.label), 'Channel count mismatch with elec.');

    a1 = d.(factA) == lvlA(1);
    a2 = d.(factA) == lvlA(2);

    if ~isInteraction
        wA = mean(c.data(:, :, a1), 3);
        wB = mean(c.data(:, :, a2), 3);
    elseif isThreeWay
        b1 = d.(factB) == lvlB(1);
        b2 = d.(factB) == lvlB(2);
        c1 = d.(factC) == lvlC(1);
        c2 = d.(factC) == lvlC(2);
        % A x B interaction (difference of differences) at each level of C
        ab = @(m) (mean(c.data(:,:,a1 & b1 & m), 3) - mean(c.data(:,:,a2 & b1 & m), 3)) ...
                - (mean(c.data(:,:,a1 & b2 & m), 3) - mean(c.data(:,:,a2 & b2 & m), 3));
        wA = ab(c1);
        wB = ab(c2);
    else
        b1 = d.(factB) == lvlB(1);
        b2 = d.(factB) == lvlB(2);
        % effect of A within each level of B (other factors averaged)
        wA = mean(c.data(:,:,a1 & b1), 3) - mean(c.data(:,:,a2 & b1), 3);
        wB = mean(c.data(:,:,a1 & b2), 3) - mean(c.data(:,:,a2 & b2), 3);
    end

    tl = struct('label', {elecData.label}, 'time', c.times, 'dimord', 'chan_time');
    dataA{i} = tl;  dataA{i}.avg = wA;
    dataB{i} = tl;  dataB{i}.avg = wB;
end

% ---- neighbours --------------------------------------------------------
cfgN = [];
cfgN.method   = 'triangulation';
cfgN.elec = elecData;
cfgN.feedback = 'no';
neighbours = ft_prepare_neighbours(cfgN);
fprintf('[3] neighbours done\n');

% ---- statistics --------------------------------------------------------
cfg = [];
cfg.channel          = 'all';
cfg.latency          = [0 max(dataA{1}.time)];
cfg.parameter        = 'avg';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 2000;
cfg.neighbours       = neighbours;
cfg.design           = [1:n 1:n; ones(1,n) 2*ones(1,n)];
cfg.uvar             = 1;
cfg.ivar             = 2;

stat = ft_timelockstatistics(cfg, dataA{:}, dataB{:});
fprintf('[3] neighbours done\n');

% ---- report ------------------------------------------------------------
fprintf('\n=== %s | %s | n = %d ===\n', label, eventName, n);
for pol = {'pos','neg'}
    f = [pol{1} 'clusters'];
    if isfield(stat, f) && ~isempty(stat.(f))
        for k = 1:numel(stat.(f))
            fprintf('  %s cluster %d: p = %.4f\n', pol{1}, k, stat.(f)(k).prob);
        end
    else
        fprintf('  no %s clusters\n', pol{1});
    end
end
   cfgGA = []; cfgGA.parameter = 'avg';
   GA_A = ft_timelockgrandaverage(cfgGA, dataA{:});
   GA_B = ft_timelockgrandaverage(cfgGA, dataB{:});
   cfgM = []; cfgM.operation = 'subtract'; cfgM.parameter = 'avg';
   GA_diff = ft_math(cfgM, GA_A, GA_B);
   save(fullfile(outFolder, sprintf('unfold_stat_%s_%s.mat', eventName, label)), 'stat', 'GA_diff');
   
save(fullfile(outFolder, sprintf('unfold_stat_%s_%s.mat', eventName, label)), 'stat');
end
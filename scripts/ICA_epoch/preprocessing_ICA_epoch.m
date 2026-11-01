%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EEG data preprocessing with Santina 
%%% Author: Santina Kemper 
%%% Last update: 05.10.2026
%%% V2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Important information before running 
% the other scripts use fieldtrip which does not go well with eeglab used
% here: so if you want to run this script again after running grand average,
% run in command window: 
% restoredefaultpath; rehash toolboxcache;
% cd('/Users/santinakemper/Documents/eeglab2026.0.0/'); eeglab; 
% which pop_loadbv (if we get a path it works, if it does not work we have to download the plugin again)
% add to path not chnage because need eeglab 

%% ========================== Structure ===================================
% Structure of this preprocessing pipeline 
% 1. Read in data and chanlocFile for electrode placement 
% 2. Filter with bandpass 0.1-30Hz 
% 3. assign electrode locations
% 4. Re-reference and remove eye channels 
% 5. ICA and automatic rejection 
% 6. Parse triggers and readjust for correct epoching
% 7. Epoching and baseline correction 
% 8. Automatic artifact rejection
% 9. Trial infor and single trial amplitudes
% 10. FieldTrip for cluster based permutation

%% ======================== set up ========================================
% to clear workspace and command window run:
%clc; clear all;

addpath(fileparts(mfilename('fullpath'))); % always find this script's own folder, regardless of cd

cd('/eeglab2026.0.0/'); % here we add the eeglab plugin 
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % initialise eeglab 

%% ======================== Configurations ================================
cfgR.rawFolder   = 'Path_rawdata'; % it seems to change a little depending on how you access the server?
cfgR.outFolder   = 'Path_preprocesseddata';

% cfgR.originalMastoidLabel = 'incase we ever actually have two'; 
% check: average mastoid channel (in Valentina's script it's 0?) 
cfgR.icArtifactThresh = 0.8; % thresholds for automatic rejection --> if more than 80% sure that this is not brain or other, reject
cfgR.icOtherThresh    = 0.1; % if other categories in other other than brain make up 10% of the source, reject

cfgR.epochWin  = [-0.1 1.0]; % window for epoching (took the numbers from Valentina's script) 
cfgR.ampThresh = [-100 100]; % uV
cfgR.jpLocalSD    = 3;   % per-channel joint-probability threshold
cfgR.jpGlobalSD   = 3;   % (computed but effectively overwritten by minJPChannels for rejection)
cfgR.kurtLocalSD  = 3;   % recorded only, not used for rejection
cfgR.kurtGlobalSD = 3;
cfgR.minJPChannels = 3;

cfgR.windows.N400 = [100 500];
cfgR.windows.Nref = [300 800];
cfgR.windows.P600 = [500 800];
cfgR.windows.ELAN = [100 300];
cfgR.roi.N400 = {'Cz','C3','C4','CP1','CP2','Pz','P3','P4'}; 
cfgR.roi.Nref = {'F3','FC1','FC5','C3'};
cfgR.roi.P600 = {'Pz','P3','P4','CP1','CP2'};
cfgR.roi.ELAN = {'Fp1','F3','FC5','FC1'};

cfgR.reviewEpochs = true;

subjectNums = 50:50; % honestly takes some time so maybe do in bunches 
% subjectNums = 13:24;
% subjectNums = 25:36;
% subjectNums = 37:48;

%% ======================== Start loop ====================================
failedSubjects = {}; % set this up to get info about processing issues (like loss of connection to the server)
for s = subjectNums
    subjectID = sprintf('P%d', s);  % creates a participant ID following the pattern VP1-VP48
    vhdr_file = sprintf('VP%d.vhdr', s);% creates the vhdr number 
    fprintf('\n === Processing %s (vhdr file: %s) === \n', subjectID, vhdr_file); % printed to give us information about where in the process we are 

    try
        %% Load + filter
        % first load the vhdr file and eeg file 
        EEG = pop_loadbv(cfgR.rawFolder, vhdr_file);
        EEG.setname = [subjectID '_raw']; % give new eeg data set the name of the participant + raw
        EEG = pop_eegfiltnew(EEG, 0.1, 30); % filter with a bandpass filter between 0.1 and 30 Hz 

        %% Channel locations
        chanlocs_full = bvef_to_chanlocs('channel locations/CMA-26-X1.bvef'); % this is the channel locations file
        EEG.chanlocs = merge_locations_by_label(EEG.chanlocs, chanlocs_full); % give EEG file information about the locations derived from the file + label 
        EEG = eeg_checkchanlocs(EEG); % check channels 

        figure; topoplot([], EEG.chanlocs, 'style', 'blank', 'electrodes', 'labels'); % plot them 

        %% Remove EOG channels, re-reference
        % remove eye channels and re-reference but keep copy
        eogChans = {'HEOGL','HEOGR','VEOGRU','VEOGRD','VEOGLU','VEOGLD'};
        eogIdx   = find(ismember({EEG.chanlocs.labels}, eogChans));
        eogData  = double(EEG.data(eogIdx,:));  
        EEG = pop_select(EEG, 'nochannel', {'HEOGL','HEOGR','VEOGRU','VEOGRD','VEOGLU','VEOGLD'});
        
        nSamples = size(EEG.data, 2); % set up matrix with channel x time sample rows and columns 
        EEG.data(end+1, :) = zeros(1, nSamples); % create a dummy channel 33 with all zeros because we reference the data to itself 
        EEG.nbchan = EEG.nbchan + 1;  % add to channel info in our eeg data 
        EEG.chanlocs(end+1).labels = 'REF_dummy'; % name it dummy, because we use it assume our online reference with A2
 
        EEG = pop_reref(EEG, {'A1', 'REF_dummy'}); % this does mathemtically this: A1-0/2 
        EEG = pop_select(EEG, 'nochannel', {'REF_dummy'}); % and then we reject the dummy reference again

        %% ICA + automatic component rejection
        % I do this automatically but it wouldn't hurt to look at the data
        % before hand 
        EEG = pop_runica(EEG, 'extended', 1); % this tells MATLAB to use the extended infomax ICA because EEG has sub-gaussian and hyper-guassian elements and it can deal with this
        EEG = pop_iclabel(EEG, 'default'); % use default labels (muscle, brain, eye etc) 
        %%
        [autoReject, icDecisions] = auto_reject_ica(EEG, cfgR.icArtifactThresh, cfgR.icOtherThresh); % this is the auto rejection set up found in the corresponding file 
        writetable(icDecisions, fullfile(cfgR.outFolder, sprintf('%s_ica_decisions.csv', subjectID))); % creates a csv file with all information about rejections and retainments 
        icAct = (EEG.icaweights * EEG.icasphere) * double(EEG.data(EEG.icachansind,:));
        C = corr(icAct', eogData');              % components x EOG channels
        maxCorr = max(abs(C), [], 2);
        eogFlag = find(maxCorr >= 0.3);          % starting threshold, check by eye
        
        rejIdx = autoReject;
        if islogical(rejIdx), rejIdx = find(rejIdx); end
        missed = setdiff(eogFlag, rejIdx);
        fprintf('%s: %d components correlate with EOG; %d of them NOT in autoReject: %s\n', ...
            subjectID, numel(eogFlag), numel(missed), mat2str(missed(:)'));

        EEG = pop_subcomp(EEG, autoReject, 0); 

        pop_saveset(EEG, 'filename', sprintf('%s_ica_processed.set', subjectID), 'filepath', cfgR.outFolder); % save file 

        %% Parse triggers -> relabel critical-trial events
        % to save the data with all relevant information we neet to relabel
        % them (function below), reminder: older versions do not support
        % functions at the buttom 
        [EEG, trialCount] = relabel_critical_trials(EEG); 
        fprintf('%s: Found %d critical trials (expected 96)\n', subjectID, trialCount); % check if I actually recorded all critical trials
        if trialCount < 90
            warning('%s: Only %d critical trials found -- check this subject''s log/recording.', ... % this indicates whether I messed up or not
                subjectID, trialCount);
        end

        %% Epoch + baseline correct
        % baseline correct and epoch based on the function below (setting
        % up the correct epochs and information about them)
        % for determiner, nouns and adjective all 100ms before start of
        % det
       
        % Apply baseline defined relative to the determiner onset
        [EEG_det1, EEG_adj1, EEG_noun1] = ...
            epoch_with_det_baseline(EEG, cfgR.epochWin);

        pop_saveset(EEG_det1,  'filename', sprintf('%s_det1_epoch.set',  subjectID), 'filepath', cfgR.outFolder); % save again so hopefully nothing is lost 
        pop_saveset(EEG_adj1,  'filename', sprintf('%s_adj1_epoch.set',  subjectID), 'filepath', cfgR.outFolder);
        pop_saveset(EEG_noun1, 'filename', sprintf('%s_noun1_epoch.set', subjectID), 'filepath', cfgR.outFolder);

        %% flag first, reject once at the end
        % amplitude threshold, joint probability, kurtosis combined 
        nBefore_det1  = EEG_det1.trials; % number of all trials before rejection
        nBefore_adj1  = EEG_adj1.trials;
        nBefore_noun1 = EEG_noun1.trials;
        
        [EEG_det1,  decisionTable_det1]  = auto_reject_epochs(EEG_det1,  cfgR);
        [EEG_adj1,  decisionTable_adj1]  = auto_reject_epochs(EEG_adj1,  cfgR);
        [EEG_noun1, decisionTable_noun1] = auto_reject_epochs(EEG_noun1, cfgR);

        writetable(decisionTable_det1,  fullfile(cfgR.outFolder, sprintf('%s_det1_rejdecisions.csv',  subjectID)));
        writetable(decisionTable_adj1,  fullfile(cfgR.outFolder, sprintf('%s_adj1_rejdecisions.csv',  subjectID)));
        writetable(decisionTable_noun1, fullfile(cfgR.outFolder, sprintf('%s_noun1_rejdecisions.csv', subjectID))); 

        pop_saveset(EEG_det1,  'filename', sprintf('%s_det1_crit.set',  subjectID), 'filepath', cfgR.outFolder); % save again so hopefully nothing is lost 
        pop_saveset(EEG_adj1,  'filename', sprintf('%s_adj1_crit.set',  subjectID), 'filepath', cfgR.outFolder);
        pop_saveset(EEG_noun1, 'filename', sprintf('%s_noun1_crit.set', subjectID), 'filepath', cfgR.outFolder);

        %% Trial info + single-trial amplitudes (for the LMM)
        % trial info contains all the information about the conditions, the
        % N400, P600, ELAN and Nref 
        TdetInfo  = computeTrialInfo(EEG_det1,  'det1',  subjectID, cfgR.windows, cfgR.roi);
        TadjInfo  = computeTrialInfo(EEG_adj1,  'adj1',  subjectID, cfgR.windows, cfgR.roi);
        TnounInfo = computeTrialInfo(EEG_noun1, 'noun1', subjectID, cfgR.windows, cfgR.roi);


        writetable(TdetInfo,  fullfile(cfgR.outFolder, sprintf('%s_det1_trialinfo.csv',  subjectID))); % creates a csv file for each particiapnt which we later use for analysis
        writetable(TadjInfo,  fullfile(cfgR.outFolder, sprintf('%s_adj1_trialinfo.csv',  subjectID)));
        writetable(TnounInfo, fullfile(cfgR.outFolder, sprintf('%s_noun1_trialinfo.csv', subjectID)));

        %% FieldTrip export (for cluster-based permutation)
        % this is important for the cluster-based permutation because we
        % want to keep all the time information 
        ftDet1  = eeg2fieldtrip_raw(EEG_det1,  TdetInfo);
        ftAdj1  = eeg2fieldtrip_raw(EEG_adj1,  TadjInfo);
        ftNoun1 = eeg2fieldtrip_raw(EEG_noun1, TnounInfo);

        save(fullfile(cfgR.outFolder, sprintf('%s_det1_ft.mat',  subjectID)), 'ftDet1'); % writes one file for each participant
        save(fullfile(cfgR.outFolder, sprintf('%s_adj1_ft.mat',  subjectID)), 'ftAdj1');
        save(fullfile(cfgR.outFolder, sprintf('%s_noun1_ft.mat', subjectID)), 'ftNoun1');

        fprintf('=== %s done ===\n', subjectID);

    catch ME
        warning('%s FAILED: %s', subjectID, ME.message); % if something goes wrong, the loop breaks and it tells me about it
        failedSubjects{end+1} = subjectID; %#ok<SAGROW>
    end
end

if ~isempty(failedSubjects)
    fprintf('\n Subjects that FAILED and need manual attention: %s\n', strjoin(failedSubjects, ', ')); % this means mostly check server connection and whether we actually connected to the folder (not just vpn)
else
    fprintf('\nAll subjects processed successfully.\n');
end

%% ========================================================================
%% Local functions
%% ========================================================================
% Disclaminer: these functions were written using the help of AI (Claude
% Sonnet-5) 

function [EEGepoched, decisionTable] = auto_reject_epochs(EEGepoched, cfgR)

% ============================================================
% AUTOMATIC EPOCH REJECTION
% Reject epoch if:
%   1. Any channel exceeds amplitude threshold
%   2. Joint probability flags >= cfgR.minJPChannels
% ============================================================

nEpochs = EEGepoched.trials;
nChan   = EEGepoched.nbchan;

% 1. Amplitude
EEGepoched = pop_eegthresh( ...
    EEGepoched, ...
    1, ...
    1:nChan, ...
    cfgR.ampThresh(1), ...
    cfgR.ampThresh(2), ...
    EEGepoched.xmin, ...
    EEGepoched.xmax, ...
    0, 0);

rejThresh = logical(EEGepoched.reject.rejthresh(:));


% 2. Joint probability
EEGepoched = pop_jointprob( ...
    EEGepoched, ...
    1, ...
    1:nChan, ...
    cfgR.jpLocalSD, ...
    cfgR.jpGlobalSD, ...
    0, 0, 0);

rejJP = logical(EEGepoched.reject.rejjp(:));


% 3. Kurtosis
EEGepoched = pop_rejkurt( ...
    EEGepoched, ...
    1, ...
    1:nChan, ...
    cfgR.kurtLocalSD, ...
    cfgR.kurtGlobalSD, ...
    0, 0, 0);

rejKurt = logical(EEGepoched.reject.rejkurt(:));


% 4. Number of JP-flagged channels
jpMatrix = false(nChan, nEpochs);

if isfield(EEGepoched.reject, 'rejjpE') && ...
        ~isempty(EEGepoched.reject.rejjpE)

    jpMatrix = logical(EEGepoched.reject.rejjpE);

end

nJPChannels = sum(jpMatrix, 1)';
nJPChannels = nJPChannels(:);


% 5. Final rejection rule
rejJPsevere = nJPChannels >= cfgR.minJPChannels;
rejJPsevere = rejJPsevere(:);

rejAll = rejThresh | rejJPsevere;
rejAll = rejAll(:);

% 6. Maximum absolute amplitude
maxAbsAmp = zeros(nEpochs, 1);

for ep = 1:nEpochs
    maxAbsAmp(ep) = max(abs(EEGepoched.data(:,:,ep)), [], 'all');
end

% 7. Reasons
reason = repmat({''}, nEpochs, 1);

for ep = 1:nEpochs

    reasons = {};

    if rejThresh(ep)
        reasons{end+1} = 'Amplitude';
    end

    if rejJPsevere(ep)
        reasons{end+1} = sprintf( ...
            'JointProb_%dChannels', nJPChannels(ep));
    end

    if isempty(reasons)
        reason{ep} = 'Keep';
    else
        reason{ep} = strjoin(reasons, ' + ');
    end

end
fprintf('\n--- DIMENSIONS BEFORE TABLE ---\n');
fprintf('nEpochs       = %d\n', nEpochs);
fprintf('Epoch         = %s\n', mat2str(size((1:nEpochs)')));
fprintf('rejThresh     = %s\n', mat2str(size(rejThresh)));
fprintf('rejJP         = %s\n', mat2str(size(rejJP)));
fprintf('rejKurt       = %s\n', mat2str(size(rejKurt)));
fprintf('nJPChannels   = %s\n', mat2str(size(nJPChannels)));
fprintf('maxAbsAmp     = %s\n', mat2str(size(maxAbsAmp)));
fprintf('rejJPsevere   = %s\n', mat2str(size(rejJPsevere)));
fprintf('rejAll        = %s\n', mat2str(size(rejAll)));
fprintf('reason        = %s\n', mat2str(size(reason)));
fprintf('--------------------------------\n');

% 8. Decision table
decisionTable = table( ...
    (1:nEpochs)', ...
    rejThresh, ...
    rejJP, ...
    rejKurt, ...
    nJPChannels, ...
    maxAbsAmp, ...
    rejAll, ...
    reason, ...
    'VariableNames', { ...
    'Epoch', ...
    'AmplitudeFlag', ...
    'JPFlag', ...
    'KurtosisFlag', ...
    'JPChannels', ...
    'MaxAbsAmplitude', ...
    'Rejected', ...
    'Reason'});

% 9. Print summary
fprintf('\n========================================\n');
fprintf('AUTOMATIC EPOCH REJECTION\n');
fprintf('========================================\n');

fprintf('Total epochs:              %d\n', nEpochs);
fprintf('Amplitude flagged:         %d\n', sum(rejThresh));
fprintf('JP flagged:                %d\n', sum(rejJP));
fprintf('JP severe:                 %d\n', sum(rejJPsevere));
fprintf('Kurtosis flagged:          %d\n', sum(rejKurt));
fprintf('Final rejected:            %d\n', sum(rejAll));
fprintf('Final retained:            %d\n', nEpochs - sum(rejAll));
fprintf('Retention:                 %.1f%%\n', ...
    100 * (1 - mean(rejAll)));

% 10. Actually remove epochs
[rejFinal, decisionTable] = review_epoch_rejections( ...
    EEGepoched, rejAll, decisionTable, cfgR);

% Remove stale ICA activation matrix before changing the number of epochs
EEGepoched.icaact = [];

% Now remove rejected epochs
EEGepoched = pop_rejepoch(EEGepoched, rejFinal, 0);

end

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

function ft = eeg2fieldtrip_raw(EEG, trialInfoTable)
    ft.label = {EEG.chanlocs.labels}';
    ft.fsample = EEG.srate;
    ft.trial = cell(1, EEG.trials);
    ft.time  = cell(1, EEG.trials);
    tvec = EEG.times / 1000;
    for e = 1:EEG.trials
        ft.trial{e} = EEG.data(:,:,e);
        ft.time{e}  = tvec;
    end
    patternNum = nan(height(trialInfoTable),1);
    patternMap = containers.Map({'A','B','C','D'}, {1,2,3,4});
    for r = 1:height(trialInfoTable)
        p = trialInfoTable.pattern{r};
        if isKey(patternMap, p), patternNum(r) = patternMap(p); end
    end
    ft.trialinfo = [trialInfoTable.item, trialInfoTable.mention, ...
        double(trialInfoTable.modified), patternNum, ...
        double(trialInfoTable.informative), trialInfoTable.block];
end

function chanlocsOut = bvef_to_chanlocs(bvefFile)

    doc = xmlread(bvefFile);
    nodes = doc.getElementsByTagName('Electrode');
    n = nodes.getLength();

    labels = cell(n,1);
    Theta = nan(n,1); Phi = nan(n,1); Radius = nan(n,1);
    valid = false(n,1);

    for i = 0:n-1
        node = nodes.item(i);
        try
            labels{i+1} = char(node.getElementsByTagName('Name').item(0).getFirstChild.getData);
            Theta(i+1)  = str2double(node.getElementsByTagName('Theta').item(0).getFirstChild.getData);
            Phi(i+1)    = str2double(node.getElementsByTagName('Phi').item(0).getFirstChild.getData);
            Radius(i+1) = str2double(node.getElementsByTagName('Radius').item(0).getFirstChild.getData);
            valid(i+1)  = true;
        catch
            valid(i+1) = false;  
        end
    end

    labels = labels(valid); Theta = Theta(valid); Phi = Phi(valid); Radius = Radius(valid);

    % --- Brainstorm's verified spherical -> Cartesian conversion ---
    phi_rad   = deg2rad(180 + Phi);
    theta_rad = deg2rad(90 - Theta);
    sx = Radius .* cos(theta_rad) .* cos(phi_rad);
    sy = Radius .* cos(theta_rad) .* sin(phi_rad);
    sz = Radius .* sin(theta_rad);

    X = -sy;        % toward nose
    Y =  sx;        % toward left ear
    Z =  sz + 0.05; % toward vertex

    chanlocsOut = struct('labels', labels(:), ...
        'X', num2cell(X(:)), ...
        'Y', num2cell(Y(:)), ...
        'Z', num2cell(Z(:)));


    chanlocsOut = convertlocs(chanlocsOut, 'cart2all');

    fprintf('bvef_to_chanlocs: parsed %d electrode positions from %s\n', ...
        numel(labels), bvefFile);

    fprintf('bvef_to_chanlocs: parsed %d electrode positions from %s\n', ...
        numel(labels), bvefFile);
end

function [EEG_det1, EEG_adj1, EEG_noun1] = ...
    epoch_with_det_baseline(EEG, epochWin)

% ============================================================
% EPOCH_WITH_DET_BASELINE
%
% Epoch DET1, ADJ1 and N1 separately.
%
% For every critical trial:
%
%   baseline = 100 ms immediately preceding DET1 onset
%
% The SAME trial-specific baseline is then subtracted from:
%
%   1. DET1 epoch
%   2. ADJ1 epoch
%   3. N1 epoch
%
% The baseline is calculated from the continuous EEG before
% epoching. This avoids using the immediately preceding
% auditory word as the baseline for ADJ1 or N1.
% ============================================================

    fs = EEG.srate;

    
    % 1. Find all DET1 critical events
    detEvents = find(strcmp({EEG.event.type}, 'DET1_crit'));

    nTrials = numel(detEvents);

    fprintf('\n');
    fprintf('DET-baseline correction\n');
    fprintf('-----------------------\n');
    fprintf('Critical trials found: %d\n', nTrials);

    % 2. Calculate one baseline per critical trial
    baselineValues = nan(EEG.nbchan, nTrials);

    baselineValid = false(1, nTrials);

    for t = 1:nTrials

        detLatency = EEG.event(detEvents(t)).latency;
        % Take the 100 ms immediately preceding DET onset.
        baseStart = round(detLatency - 0.100 * fs);
        baseEnd   = round(detLatency) - 1;
        if baseStart < 1 || baseEnd > EEG.pnts

            warning(['Critical trial %d: DET baseline is outside ', ...
                     'the continuous EEG.'], t);

            continue;
        end

        % Mean across time, separately for every channel
        baselineValues(:, t) = mean( ...
            EEG.data(:, baseStart:baseEnd), ...
            2, 'omitnan');

        baselineValid(t) = true;

    end

    fprintf('Valid DET baselines: %d/%d\n', ...
        sum(baselineValid), nTrials);

    % 3. Epoch WITHOUT baseline correction
    EEG_det1 = pop_epoch( ...
        EEG, {'DET1_crit'}, epochWin);

    EEG_adj1 = pop_epoch( ...
        EEG, {'ADJ1_crit'}, epochWin);

    EEG_noun1 = pop_epoch( ...
        EEG, {'N1_crit'}, epochWin);

    % 4. Apply the corresponding DET baseline
    EEG_det1 = apply_det_baseline( ...
        EEG_det1, baselineValues, baselineValid);

    EEG_adj1 = apply_det_baseline( ...
        EEG_adj1, baselineValues, baselineValid);

    EEG_noun1 = apply_det_baseline( ...
        EEG_noun1, baselineValues, baselineValid);

    fprintf('DET1, ADJ1 and N1 baseline correction complete.\n\n');

end


function EEG_out = apply_det_baseline( ...
    EEG_in, baselineValues, baselineValid)

    EEG_out = EEG_in;

    nEpochs = EEG_in.trials;

    for e = 1:nEpochs

        if ~isfield(EEG_in.epoch, 'eventcriticalTrialID')

            error(['Epoch structure does not contain ', ...
                   'eventcriticalTrialID. Check that ', ...
                   'relabel_critical_trials was run before ', ...
                   'epoching.']);

        end

        trialID = EEG_in.epoch(e).eventcriticalTrialID;

        
        % Convert possible cell/array representation to scalar
        if iscell(trialID)

            % Keep non-empty values
            trialID = trialID(~cellfun(@isempty, trialID));

            if isempty(trialID)
                warning('Epoch %d has no criticalTrialID. Skipping.', e);
                continue;
            end

            trialID = trialID{1};

        end

        if numel(trialID) > 1
            trialID = trialID(1);
        end

        trialID = double(trialID);

        % Validate trial ID
        if isnan(trialID) || ...
           trialID < 1 || ...
           trialID > size(baselineValues, 2)

            warning(['Epoch %d has invalid criticalTrialID ', ...
                     '(%s). Skipping.'], ...
                     e, mat2str(trialID));

            continue;
        end

        trialID = round(trialID);

       if ~baselineValid(trialID)

            warning(['Epoch %d: baseline for critical trial ', ...
                     '%d was invalid. Skipping.'], ...
                     e, trialID);

            continue;
        end

        % Subtract trial-specific DET baseline
        
        EEG_out.data(:,:,e) = ...
            EEG_out.data(:,:,e) - baselineValues(:,trialID);

    end

    EEG_out = eeg_checkset(EEG_out);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EEG data preprocessing: Deconvolution correction 
%%% Author: Santina Kemper 
%%% Last update: 01.10.2026
%%% V1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================== Structure ===================================
% opens unfold and runs deconvolution 
% this script was written with the help of Claude Sonnet 5.5
%% ========================== Set up ===================================

addpath(fileparts(mfilename('fullpath'))); % always find this script's own folder, regardless of cd
cd('/eeglab2026.0.0/'); % here we add the eeglab plugin 
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; 
run('/unfold/init_unfold.m');

cfgR.inFolder   = 'ica_preprocessed_data_path'; 
cfgR.outFolder   = 'unfold_preprocessed_path'; 
subjectNums = 1:17;
cfgR.epochWin  = [-0.1 1.0]; % window for epoching (took the numbers from Valentina's script) 

cfgR.windows.N400 = [100 500];
cfgR.windows.Nref = [300 800];
cfgR.windows.P600 = [500 800];
cfgR.windows.ELAN = [100 300];
cfgR.roi.N400 = {'Cz','C3','C4','CP1','CP2','Pz','P3','P4'}; 
cfgR.roi.Nref = {'F3','FC1','FC5','C3'};
cfgR.roi.P600 = {'Pz','P3','P4','CP1','CP2'};
cfgR.roi.ELAN = {'Fp1','F3','FC5','FC1'};

%% ========================== Deconvolution ===============================
for s = subjectNums
    subjectID = sprintf('P%d', s);
    EEG = pop_loadset('filename', sprintf('%s_ica_processed.set', subjectID), 'filepath', cfgR.inFolder);

    [EEG, trialCount] = relabel_critical_trials(EEG); 
        fprintf('%s: Found %d critical trials (expected 96)\n', subjectID, trialCount); % check if I actually recorded all critical trials
        if trialCount < 90
            warning('%s: Only %d critical trials found -- check this subject''s log/recording.', ... % this indicates whether I messed up or not
                subjectID, trialCount);
        end

    % --- design: one formula per event type ---
    % det1 and noun1: all three predictors vary and can interact.
    % adj1: 'modified' is constant (always 1) -- including it would be
    % collinear/degenerate, so it's dropped from that event's formula
    cfgDesign = [];
    cfgDesign.eventtypes = {'DET1_crit', 'ADJ1_crit', 'N1_crit'};
    cfgDesign.formula = { ...
        'y ~ 1 + cat(mention)*cat(modified)*cat(informative)', ...
        'y ~ 1 + cat(mention)*cat(informative)', ...
        'y ~ 1 + cat(mention)*cat(modified)*cat(informative)'};
    EEG = uf_designmat(EEG, cfgDesign);

    % --- time-expand: same window per event as in existing epoching
    cfgTimeshift = [];
    cfgTimeshift.timelimits = cfgR.epochWin; 
    EEG = uf_timeexpandDesignmat(EEG, cfgTimeshift);

    % --- continuous artifact handling (replaces epoch rejection) ---
    winrej = uf_continuousArtifactDetect(EEG, 'amplitudeThreshold', 150);

    badSamples = false(1, EEG.pnts);
    
    for k = 1:size(winrej, 1)
        firstSample = max(1, round(winrej(k, 1)));
        lastSample  = min(EEG.pnts, round(winrej(k, 2)));
    
        if firstSample <= lastSample
            badSamples(firstSample:lastSample) = true;
        end
    end
    
    pctBad = 100 * mean(badSamples);
    
    fprintf('%s: %d marked segments, %.1f%% of samples flagged\n', ...
        subjectID, size(winrej, 1), pctBad);

    EEG = uf_continuousArtifactExclude(EEG, struct('winrej', winrej));

    % --- fit the deconvolution model (continuous data, not epoched) ---
    EEG = uf_glmfit(EEG);

    % --- export singel trial data ---
    chanIdx = find(ismember({EEG.chanlocs.labels}, ...
    unique([cfgR.roi.N400, cfgR.roi.Nref, cfgR.roi.P600, cfgR.roi.ELAN])));
    Tst = export_single_trials(EEG, subjectID, badSamples, cfgR.epochWin, chanIdx, 20);
    writetable(Tst, fullfile(cfgR.outFolder, sprintf('%s_single_trials.csv', subjectID)));
    
    % --- condense into per-condition modeled waveforms ---
    ufresult = uf_condense(EEG);

    save(fullfile(cfgR.outFolder, sprintf('%s_unfold_result.mat', subjectID)), 'ufresult');

    save_subject_cells(ufresult, subjectID, cfgR.outFolder, [-0.1 0]);
    
end


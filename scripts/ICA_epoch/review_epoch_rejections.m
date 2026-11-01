function [rejFinal, decisionTable] = review_epoch_rejections( ...
    EEGepoched, rejAll, decisionTable, cfgR)

% REVIEW_EPOCH_REJECTIONS
%
% 1. Show automatically flagged epochs in EEGLAB.
% 2. Close the EEGLAB plot after visual inspection.
% 3. Enter epochs to KEEP despite automatic rejection.
% 4. Enter epochs to REJECT despite automatic retention.
% 5. Return the final rejection vector.
%
% No epochs are removed here.
% The calling function performs pop_rejepoch() afterwards.

    % =============================================================
    % 0. Skip review if disabled
    % =============================================================

    if ~isfield(cfgR, 'reviewEpochs') || ~cfgR.reviewEpochs
        rejFinal = logical(rejAll(:));
        return
    end

    % =============================================================
    % 1. Basic dimensions
    % =============================================================

    nEpochs = EEGepoched.trials;
    nChan   = EEGepoched.nbchan;

    rejAll = logical(rejAll(:));

    if numel(rejAll) ~= nEpochs
        error(['review_epoch_rejections: rejAll contains %d values, ' ...
               'but EEGepoched contains %d epochs.'], ...
               numel(rejAll), nEpochs);
    end

    if height(decisionTable) ~= nEpochs
        error(['review_epoch_rejections: decisionTable has %d rows, ' ...
               'but EEGepoched contains %d epochs.'], ...
               height(decisionTable), nEpochs);
    end

    % =============================================================
    % 2. Give EEGLAB the automatic epoch decisions for display
    % =============================================================

    EEGepoched.reject.rejmanual = rejAll;

    % EEGLAB expects channel x epoch information here
    EEGepoched.reject.rejmanualE = repmat(rejAll', nChan, 1);

    % =============================================================
    % 3. Tell user what is happening
    % =============================================================

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('%s | MANUAL EPOCH REVIEW\n', EEGepoched.setname);
    fprintf('============================================================\n');
    fprintf('Total epochs:              %d\n', nEpochs);
    fprintf('Automatically flagged:     %d\n', sum(rejAll));
    fprintf('Automatically retained:    %d\n', nEpochs - sum(rejAll));
    fprintf('\n');
    fprintf('Inspect the epochs in the EEGLAB window.\n');
    fprintf('The automatically flagged epochs are highlighted.\n');
    fprintf('\n');
    fprintf('When you are finished inspecting them, CLOSE the EEGLAB window.\n');
    fprintf('The Command Window will then ask for your manual decisions.\n');
    fprintf('============================================================\n');

    % =============================================================
    % 4. Open EEGLAB plot
    % =============================================================

    pop_eegplot(EEGepoched, 1, 1, 0);

    % =============================================================
    % 5. IMPORTANT:
    % These variables MUST be created after the plot closes.
    % =============================================================

    keepOverride = input( ...
        ['Epoch #s to KEEP despite the automatic flag, ' ...
         'e.g. [3 17] (Enter = none): '] );

    rejOverride = input( ...
        ['Epoch #s to REJECT that were NOT flagged, ' ...
         'e.g. [5 22] (Enter = none): '] );

    % =============================================================
    % 6. Start with automatic decision
    % =============================================================

    rejFinal = rejAll;

    % =============================================================
    % 7. Manual KEEP
    % =============================================================

    if ~isempty(keepOverride)

        keepOverride = unique(keepOverride(:)');

        if any(keepOverride < 1 | ...
               keepOverride > nEpochs | ...
               keepOverride ~= round(keepOverride))

            error('Invalid epoch number in keepOverride.');
        end

        rejFinal(keepOverride) = false;

        decisionTable.Reason(keepOverride) = ...
            strcat(decisionTable.Reason(keepOverride), ...
                   ' +ManualKeep');

    end

    % =============================================================
    % 8. Manual REJECT
    % =============================================================

    if ~isempty(rejOverride)

        rejOverride = unique(rejOverride(:)');

        if any(rejOverride < 1 | ...
               rejOverride > nEpochs | ...
               rejOverride ~= round(rejOverride))

            error('Invalid epoch number in rejOverride.');
        end

        rejFinal(rejOverride) = true;

        decisionTable.Reason(rejOverride) = ...
            strcat(decisionTable.Reason(rejOverride), ...
                   ' +ManualReject');

    end

    % =============================================================
    % 9. Final decision
    % =============================================================

    rejFinal = logical(rejFinal(:));

    decisionTable.Rejected = rejFinal;

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('FINAL EPOCH DECISION\n');
    fprintf('============================================================\n');
    fprintf('Automatically rejected:    %d\n', sum(rejAll));
    fprintf('Manually retained:         %d\n', ...
        sum(rejAll & ~rejFinal));
    fprintf('Manually rejected:         %d\n', ...
        sum(~rejAll & rejFinal));
    fprintf('FINAL rejected:            %d\n', sum(rejFinal));
    fprintf('FINAL retained:            %d\n', nEpochs - sum(rejFinal));
    fprintf('============================================================\n');

end
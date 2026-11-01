%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% EEG data preprocessing with Santina - Auto reject function 
%%% Author: Santina Kemper 
%%% Last update: 07.09.2026
%%% V2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Goal of this function 
% if OTHER in ICA is high, we need to check where the EEG could come from 
% if in OTHER Brain is a runner up, we keep Other 
% but if Brain is below 10% and another category is above 10% we reject
% other: Brain: 2%, eye 12% and Other 80% --> Reject 
% the actual thresholds can be set in the main script 
% this script was written with the assistence of Claude Sonnet 5.5

function [autoReject, decisionTable] = auto_reject_ica(EEG_ICA, artifactThresh, otherThresh)

if nargin < 2, artifactThresh = 0.8; end
if nargin < 3, otherThresh = 0.1; end

classifications = EEG_ICA.etc.ic_classification.ICLabel.classifications;
nIC = size(classifications, 1);

BRAIN = 1; MUSCLE = 2; EYE = 3; HEART = 4; LINE = 5; CHANNEL = 6; OTHER = 7;
artifactCols = [MUSCLE EYE HEART LINE CHANNEL];
labels = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};

[maxProb, maxClass] = max(classifications, [], 2);
rejectIC = false(nIC, 1);
reason = cell(nIC, 1);

for ic = 1:nIC
    probs = classifications(ic, :);

    if maxClass(ic) ~= BRAIN && maxClass(ic) ~= OTHER && maxProb(ic) > artifactThresh
        rejectIC(ic) = true;
        reason{ic} = sprintf('%s > %.0f%%', labels{maxClass(ic)}, artifactThresh*100);

    
elseif maxClass(ic) == OTHER % if other is the largest class

    brainProb = probs(BRAIN); % brain probability 
    artifactProbs = probs(artifactCols); % other probabilities 

    if brainProb >= 0.10 % above 10% brain  
        rejectIC(ic) = false; % no rejection 
        reason{ic} = sprintf( ...
            'OTHER, Brain is %.0f%% >= 10%%, therefore OTHER is kept', ...
            brainProb * 100); % give us reason for why we keep it to do sanity checks 

    elseif any(artifactProbs > otherThresh) % rejection else for other situation 
        rejectIC(ic) = true; 

        [~, worst] = max(artifactProbs);
        worstLabel = labels{artifactCols(worst)};

        reason{ic} = sprintf( ...
            'OTHER, Brain %.0f%% < 10%% and %s > %.0f%%, therefore OTHER was rejected', ...
            brainProb * 100, worstLabel, otherThresh * 100);

    else % keep if no other ICA component is large but brain also not too large 
        rejectIC(ic) = false;
        reason{ic} = sprintf( ...
            'OTHER, Brain %.0f%% < 10%%, but no artifact > %.0f%% was found, therefore OTHER was kept', ...
            brainProb * 100, otherThresh * 100);
    end
    else
        reason{ic} = sprintf('%s -> kept', labels{maxClass(ic)});
    end
end

autoReject = find(rejectIC);

decisionTable = table((1:nIC)', labels(maxClass)', maxProb, rejectIC, reason, ...
    'VariableNames', {'IC','TopClass','TopProb','Rejected','Reason'});

fprintf('Auto-rejected %d / %d components.\n', numel(autoReject), nIC);
otherRows = strcmp(decisionTable.TopClass, 'Other');
fprintf('  Of %d "Other"-labeled components, %d were rejected under the co-activation rule.\n', ...
    sum(otherRows), sum(otherRows & decisionTable.Rejected));

end

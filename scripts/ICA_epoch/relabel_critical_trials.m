function [EEG, trialCount] = relabel_critical_trials(EEG)

    condMap = containers.Map( ...
        {'109','110','111','112','113','114','115','116'}, ...
        {struct('pattern','A','informative',true),  ...
         struct('pattern','A','informative',false), ...
         struct('pattern','B','informative',true),  ...
         struct('pattern','B','informative',false), ...
         struct('pattern','C','informative',true),  ...
         struct('pattern','C','informative',false), ...
         struct('pattern','D','informative',true),  ...
         struct('pattern','D','informative',false)});

    events = strtrim(erase({EEG.event.type}', 'S'));

    nEv = numel(events);

    curBlock = NaN;
    curCondCode = '';
    curPattern = '';
    curInformative = NaN;

    trialCount = 0;

    i = 1;

    while i <= nEv

        code = events{i};

        if ismember(code, {'105','106','107','108'})
            curBlock = str2double(code) - 104;

        elseif isKey(condMap, code)

            curCondCode = code; %#ok<NASGU>

            c = condMap(code);

            curPattern = c.pattern;
            curInformative = c.informative;

        elseif strcmp(code, '200')

            itemCode = events{i+2};
            typeCode = events{i+3};

            if strcmp(typeCode, '97')
                mentionCode = events{i+4};
                wordStart = i + 5;
            else
                mentionCode = '';
                wordStart = i + 4;
            end

            wordEnd = wordStart;

            while wordEnd <= nEv && ~strcmp(events{wordEnd}, '104')
                wordEnd = wordEnd + 1;
            end

            trialIdx = wordStart:wordEnd;
            trialEvents = events(trialIdx);

            isModified = any(strcmp(trialEvents, '202'));

            if strcmp(typeCode, '97')

                % Unique ID for this critical trial
                trialCount = trialCount + 1;
                criticalTrialID = trialCount;

                for k = trialIdx

                    switch events{k}

                        case '201'
                            EEG.event(k).type = 'DET1_crit';

                        case '202'
                            EEG.event(k).type = 'ADJ1_crit';

                        case '203'
                            EEG.event(k).type = 'N1_crit';

                        otherwise
                            continue
                    end

                    EEG.event(k).item         = str2double(itemCode);
                    EEG.event(k).mention      = str2double(mentionCode);
                    EEG.event(k).modified     = isModified;
                    EEG.event(k).pattern      = curPattern;
                    EEG.event(k).informative  = curInformative;
                    EEG.event(k).block        = curBlock;

                    % Same unique ID for DET1, ADJ1 and N1 belonging
                    EEG.event(k).criticalTrialID = criticalTrialID;

                end
            end

            i = wordEnd;

        end

        i = i + 1;

    end

end
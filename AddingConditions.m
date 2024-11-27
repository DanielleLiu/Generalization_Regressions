% ------- Deprecated, use new implementation SeparateConditions.m from labtools ----------------
%Chagingparams, create new conditions in 1 trial, currently only need to
%use it fors the first block for preintervention trials.
function in= AddingConditions(in, oldConditionName, newConditions, splitBeginOrEnd, newDecription)
    warning('This function is deprecated. Please update your code to use the SeparateConditions from labTools instead.')
    if nargin < 4 || isempty(splitBeginOrEnd)
        splitBeginOrEnd = true; %default split beginning
    end

    if nargin < 5 || isempty(newDecription)
        newDecription = 'Beginning of adaptation';
    end
  
%     in.plotAvgTimeCourse(in,{'netContributionNorm2','singleStanceSpeedFastAbsANK','singleStanceSpeedSlowAbsANK'})
%     title('Before Changing Conditions');
    
    % find the index for the label singleStanceSpeedFastAbsANK
    idxfast=compareListsNested({'singleStanceSpeedFastAbsANK'},in.data.labels)==1;
    idxslow=compareListsNested({'singleStanceSpeedSlowAbsANK'},in.data.labels)==1;
    
    trialNum = in.metaData.trialsInCondition{strcmp(in.metaData.conditionName,oldConditionName)};
    columnIdxForTrialNum=find(compareListsNested({'Trial'},in.data.labels));

    fast=in.data.Data(:,idxfast);
    slow=in.data.Data(:,idxslow);
    difference=fast-slow;

    % find the index with speed difference or speed at 1 but also within the
    % current condition of interes
    idxSplit=find(abs(difference)>100 & in.data.Data(:,columnIdxForTrialNum) == trialNum);
    diffMask = (abs(difference)>100 & in.data.Data(:,columnIdxForTrialNum) == trialNum)'; %find all cases where diff in speed > 200
    thresholdStrides = 4; %if happens continuously for 4 strides: speed changed instead of noise in speed
    thresholdFramesMask = ones(1, thresholdStrides);
    % pad 0 to diffMask in case started difference at frame 1, the 2nd argument is the pattern
    % to match, find the index where the previous frame didn't have speed diff > 200, and the 
    % next 150 frames have speed diff > 200
    if splitBeginOrEnd
%         idxSplit = idxSplit(1);
        idxSplit = strfind([0,diffMask],[0,thresholdFramesMask]); %get the first stride that speeds changed
%         idxSplit = idxSplit - thresholdStrides; 
    else
%         idxSplit = idxSplit(end)+1;
        idxSplit = strfind([diffMask,0],[thresholdFramesMask,0]); %get the first stride where 4 stride later speed changed
        idxSplit = idxSplit + thresholdStrides; %shift to the first stride where speed is different now.
    end
    
    currConditionLength = length(in.metaData.conditionName);
    
    %for now handles 1 new condition at a time
    % ismember finds A in B 
%     newConditionTrialCount = 1;
    [condExist, loc] = ismember(newConditions, in.metaData.conditionName);
    condExistIdx = find(condExist);
    loc = loc(condExistIdx);

    if condExistIdx
%         Approach2: shift everything by 1
        %append an additional trial to the new condition
        newTrialsInConds = {};
        newTrialsInConds{1} = [trialNum+1 in.metaData.trialsInCondition{loc}+1];
        %for the ones following, increment by 1
        for i = loc+1:currConditionLength
            if ~isempty(in.metaData.trialsInCondition{i})
                newTrialsInConds{i-loc+1} = in.metaData.trialsInCondition{i}+1;
            else
                newTrialsInConds{i-loc+1} = [];
            end
        end
        in.metaData.trialsInCondition(loc:end) = newTrialsInConds;
        in.data = in.data.setTrialTypes({in.data.trialTypes{1:trialNum} in.data.trialTypes{trialNum+1} in.data.trialTypes{trialNum+1:end}});
        in.data.Data(idxSplit:end,columnIdxForTrialNum) = in.data.Data(idxSplit:end,columnIdxForTrialNum) +1;
    else
        [condExist, loc] = ismember(oldConditionName, in.metaData.conditionName);
        loc = loc(find(condExist));
        newTrialsInConds{1} = [trialNum+1];
        %for the ones following, increment by 1
        maxTrial = trialNum+1;
        for i = loc+1:currConditionLength
            if ~isempty(in.metaData.trialsInCondition{i})
                newTrialsInConds{i-loc+1} = in.metaData.trialsInCondition{i}+1;
                maxTrial = newTrialsInConds{i-loc+1};
                maxTrial = maxTrial(end);
            else
                newTrialsInConds{i-loc+1} = [];
            end
        end
        in.metaData.trialsInCondition = [in.metaData.trialsInCondition{1:loc}, newTrialsInConds];
        in.metaData.conditionName = [in.metaData.conditionName{1:loc}, {newConditions}, in.metaData.conditionName{loc+1:end}];
        in.metaData.conditionDescription = [in.metaData.conditionDescription{1:loc}, {newDecription}, in.metaData.conditionDescription{loc+1:end}];
        in.metaData.Ntrials = maxTrial;
        in.data = in.data.setTrialTypes({in.data.trialTypes{1:trialNum} in.data.trialTypes{trialNum} in.data.trialTypes{trialNum+1:end}});
        in.data.Data(idxSplit:end,columnIdxForTrialNum) = in.data.Data(idxSplit:end,columnIdxForTrialNum) +1;
        
%         in.metaData.conditionName{currConditionLength+1}=newConditions;
%         in.metaData.trialsInCondition{currConditionLength+1}=currMaxTrial+1;
%         in.metaData.conditionDescription{currConditionLength+1}= newDecription;
%         in.data.trialTypes{currConditionLength+1}=newTrialType;
%         in.data.Data(idxSplit,columnIdxForTrialNum)=currMaxTrial+1;
    end
    
    in.plotAvgTimeCourse(in,{'netContributionNorm2','singleStanceSpeedFastAbsANK','singleStanceSpeedSlowAbsANK'})
    title('After Adding Conditions')
end 

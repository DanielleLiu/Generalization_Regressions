function [RemovedData]=RemoveBadMuscles(normalizedGroupData,badSubjID,badMuscles )
%This is a function that change the values of the "bad" muscle to NaN. The subjects and muscle
%are hard code. Make sure you update this to your needs and that you
%removed muscle bilaterally

subjectsToPlot=[];  subjectsToPlotID=[];
    subjectsToPlot{end+1} = normalizedGroupData; %
%     subjectsToPlotID{end+1} = groupID;% from SLcode

% Bad muscles for group plots
% % %CTR
% if contains(groupID,'PATS')
%     badSubjID = {'PATS05','PATS03'}; %badSubj and muscle are index matched, if want to remove group, put group ID here
%     badMuscles = {{'sRFs', 'fRFs'},{'sRFs', 'fRFs','sHIPs','fHIPs','sVLs', 'fVLs','sVMs', 'fVMs'}}; %labels in group ID will be removed for all regression and AE computations;
% elseif contains(groupID,'PATR')
%     
% end


for idxToRemove = 1:numel(badSubjID)
    
    subjIdx = find(contains(subjectsToPlot{end}.ID, badSubjID{idxToRemove}));

        
    
    if ~isempty(subjIdx)
        
         badSubj = subjectsToPlot{end}.adaptData{subjIdx};
        
        for i = 1:numel(badMuscles{idxToRemove})
            
            
            badDataIdx=find(contains(badSubj.data.labels, {[badMuscles{idxToRemove}{i},' ']}));
            if length(badDataIdx)<12 %each muscle should have 12 sub intervals
                badDataIdxlast=badDataIdx(end)+[1:3]; %if not add 3 more intervals (rows of data) to remove. This seems rather arbitrary.
                badDataIdx= [badDataIdx; badDataIdxlast'];
            end


            badSubj.data.Data(:,badDataIdx) = nan;
            
           disp(['Removing (Setting NaN) of ' badMuscles{idxToRemove}{i} ' from Subject: ' badSubj.subData.ID])
            
        end
         subjectsToPlot{end}.adaptData{subjIdx} = badSubj;
    end
    
    
        
    
end

RemovedData=subjectsToPlot{1}; % from SL code
end

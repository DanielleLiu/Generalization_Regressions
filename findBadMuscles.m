% This script use a dialog to prompt users to enter bad muscles to remove.
% Used for quality check before running the full analysis.
%% prep data loading & set up parameters
% musclesToRemove = {{},{}};
close all; clc; clearvars -except musclesToRemove
visitNum = 'V04';
[dataPath, ~, ~, ~, ~] = setupDataPath('ParamsForGroupPlot', visitNum, '', '');
cd(dataPath) %go to where the list of params directory for data loading.

%%
subjectIDNums= [1:7 9:13 15:17 19:25];
for subIdx = length(subjectIDNums):length(subjectIDNums)
    groupID = ['AUF' num2str(subjectIDNums(subIdx),'%02d')];
    groupID
    [~, ~, resDir, subjectID, ~] = setupDataPath(groupID, visitNum, '', 'EMG');
    saveResAndFigures = false;

    %% get C, check data, add data to remove.
    removeBias = false; %to visualize data, don't remove bias to see all epochs.
    removeBadmuscles = true;
    loadDataByStep = true; %load data fresh
    gettingCs_dynamicsFeedback
    
    resDir
    if ~isfolder(resDir) %if doesn't exist, make it.
        mkdir(resDir)
    end
    saveas(fh, [resDir groupID '_checker_keyEpochs.png'])
    %% identify and populate bad muscles list.
    answer = inputdlg("Enter muscle names to remove (separated by comma, e.g.: fVLs, sVLs; Leave empty if no bad muscles.");
    if ~isempty(answer{1})
        if ~exist('musclesToRemove', 'var')
            musclesToRemove = {{},{}}; %first time initialize the variable.
        end
        musclesToRemove{1}{end+1} = [groupID visitNum];
        musclesToRemove{2}{end+1} = strtrim(split(answer{1},','))';
        % run analysis again to see the removed data.
        gettingCs_dynamicsFeedback
        saveas(fh, [resDir groupID '_checker_keyEpochs_removed' answer{1} '.png'])
    end
end

%% save the musclesToRemove
% save('musclesToRemove.mat','musclesToRemove')

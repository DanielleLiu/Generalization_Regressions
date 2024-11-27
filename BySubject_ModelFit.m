%% prep data loading & set up parameters
% musclesToRemove = {{},{}};
close all; clc; clearvars -except dataBySubjects musclesToRemove
[dataPath, ~, ~, ~, ~] = setupDataPath('ParamsForGroupPlot', 'V02', '', '');
cd(dataPath) %go to where the list of params directory for data loading.

%%
subjectIDNums= [1:7 9:13 15:17 19:25];
for subIdx = 18:length(subjectIDNums)
    close all; clc; clearvars -except dataBySubjects musclesToRemove subjectIDNums subIdx
    groupID = ['AUF' num2str(subjectIDNums(subIdx),'%02d')];
    groupID
    [~, ~, resDir, subjectID, ~] = setupDataPath(groupID, 'V02', '', 'EMG');
    saveResAndFigures = false;

    %TODO: remove muscle non symmetrical, and turn off plotting maybe to
    %get the weights.
    %% now with bad muscle list updated. Get data and run model again.
    loadDataByStep = true; %load data fresh
    removeBadmuscles = true;
    preProcessingLinearModel
    loadDataByStep = false; %avoid more data loading.
    newLabelPrefix = newLabelPrefix(end:-1:1) %reverse label orders
    removeBias = true; %get bias removed Cs.
    gettingCs_dynamicsFeedback
    model_fit_individual_muscles
    
    saveas(figure(2),[resDir, 'modelFit.png'])
    saveas(figure(3),[resDir, 'earlyPostWeights.png'])

    %% save data
    dataBySubjects.C_indv(:,:,subIdx) = C_indv;
    dataBySubjects.Ymuscles(:,:,:,subIdx) = Ymuscles;
    dataBySubjects.weights_tranp(:,:,:,subIdx) = weights_tranp;
    dataBySubjects.R2(:,subIdx) = R2'; %stride x subjects
end

%% save data
% save('IndivSub_indivMuslce_modelFitResults.mat','dataBySubjects')
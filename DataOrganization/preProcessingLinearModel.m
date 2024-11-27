%%
clear; close all; clc;
%configure parameters
adapt = false;
visitNum = 'V02';
groupID ='BATS'; %Group of interest 
if strcmpi(groupID,'AUF')
    [dataPath, ~, ~, ~, ~] = setupDataPath('ParamsForGroupPlot', visitNum, '', '');
    cd(dataPath) %go to where the list of params directory for data loading.
    load('MusclesToRemove.mat')
end
[group2, newLabelPrefix,n,subID]=creatingGroupdataWnormalizedEMG(groupID,1); % Creating the groupData normalized
%% Removing bad muscles 
%This script make sure that we always remove the same muscle for the
%different analysis 
removeBadmuscles=1;
if removeBadmuscles==1
    if exist('musclesToRemove','var')
        group2= RemovingBadMuscleToSubj(group2,musclesToRemove{1},musclesToRemove{2});
    else
        group2= RemovingBadMuscleToSubj(group2);
    end
end
%% Define epochs
strides=[-40 300]; %Number per strides per condition
if contains(groupID,{'BAT'})
    cond={'OG base','Post 1'}; %Conditions for this group
elseif contains(groupID,{'AUF'})
    if adapt
        strides=[-40 900];
        cond={'TMBase','Adaptation'}; %Conditions for this group
    else
        cond={'OGBase','OGPost'}; %Conditions for this group
    end
end
exemptFirst=[1]; %ignore inital strides
exemptLast=[5]; %Strides needed

ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Post1'}); %epochs
%% Pick muscles that you wanted to get the data from 
%%%% load and prep data
% muscleOrder={'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT'};
muscleOrder={'TA','PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
newLabelPrefix2 = defineMuscleList(muscleOrder); %List of muscle
newLabelPrefix2 = regexprep(newLabelPrefix2,'_s','s'); %removing the underscore "_"

for m=1:length(newLabelPrefix2)

    idx(m)=find(strcmp(newLabelPrefix,newLabelPrefix2(m))); %getting the index of the muscles

end

wanted_Muscles= newLabelPrefix(sort(idx)); %It needs to be the muscle form the slow leg first
newLabelPrefix= wanted_Muscles;
%% get data:
padWithNaNFlag=true;
[dataEMG,labels,allDataEMG2]=group2.getPrefixedEpochData(newLabelPrefix,ep,padWithNaNFlag); 

%Flipping EMG:
for i=1:length(allDataEMG2)
    aux=reshape(allDataEMG2{i},size(allDataEMG2{i},1),size(labels,1),size(labels,2),size(allDataEMG2{i},3));
    allDataEMG2{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
end

[~,~,dataContribs]=group2.getEpochData(ep,{'netContributionNorm2'},padWithNaNFlag);


%% Getting the regressors values
%% 
if contains(groupID,{'BAT','BAT'})
    ep=defineRegressorsDynamicsFeedback('nanmean'); 
    epochOfInterest={'Ramp','PosShort_{early}','Adaptation_{early}','Optimal','NegShort_{early}','TM base'};
elseif contains(groupID,{'AUF'})
    if contains(group2.ID{1},'V03')
        intervention = true;
    else
        intervention = false;
    end
    ep = defineEpochNirs(intervention, 'nanmean'); %thie nan is used to summarize strides in an epoch
    epochOfInterest = {'OGBase','TMBase','Adaptation_{SS}','NegShort_{la}','OGPost_{Early}','TMPost_{Early}','PosShort_{la}'};
end

padWithNaNFlag=true; %If no enough steps fill with nan, let this on

for l=1:length(epochOfInterest)
   
    ep2=defineReferenceEpoch(epochOfInterest{l},ep);
    
    [dataEMG,labels,allDataEMG]=group2.getPrefixedEpochData(newLabelPrefix,ep2,padWithNaNFlag); %Getting the data
    
    %Flipping EMG:0
    for i=1:length(allDataEMG)
        aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
        allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
        
    end
    
  regressors{l}=nanmean(allDataEMG{:},1);


end

%% Reorganize data
EMGdata=cell2mat(allDataEMG2);
regressors=cell2mat(regressors');

%%
%%Save to hdf5 format for sharing with non-Matlab users
if adapt
    name=['dynamicsData_',groupID,'_subj_', num2str(n),'_RemoveBadMuscles_', num2str(removeBadmuscles),'_Adapt_TMBase','.h5'];
else
    name=['dynamicsData_',groupID,'_subj_', num2str(n),'_RemoveBadMuscles_', num2str(removeBadmuscles),'_OGPost_OGBase','.h5'];
end

h5create(name,'/EMGdata',size(EMGdata))
h5write(name,'/EMGdata',EMGdata)

h5create(name,'/Regressors',size(regressors))
h5write(name,'/Regressors',regressors)

hdf5write(name,'/SubID',subID(:),'WriteMode','append')

hdf5write(name,'/Epochs',epochOfInterest(:),'WriteMode','append')

SLA=squeeze(cell2mat(dataContribs));
h5create(name,'/SLA',size(SLA))
h5write(name,'/SLA',SLA)
% speedDiff=[zeros(1,abs(strides(1))),ones(1,strides(2)),zeros(1,(strides(3)))];
speedDiff=[zeros(1,abs(strides(1))),ones(1,strides(2))];
h5create(name,'/speedDiff',size(speedDiff))
h5write(name,'/speedDiff',speedDiff)
breaks=[zeros(1,length(speedDiff))];
h5create(name,'/breaks',size(breaks))
h5write(name,'/breaks',breaks)
hdf5write(name,'/labels',newLabelPrefix(:),'WriteMode','append')

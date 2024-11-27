%getting Cs, this requires the current directory to be in the folder that
%contains all the params files for the relevant group.
%% Load data and Plot checkerboard for all conditions.
if exist('loadDataByStep','var') && loadDataByStep %this should only be true if running scripts step by step; 
    %default false (e.g., if the var doesn't exist, dont' load data)
    %and this variable is set somewhere outside of this script.
%     clear; close all; clc;
%     groupID ='AUF';
    [normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID);
    %the label order here is different from the previous step:
    %preProcessingLinearModel, previous step needs inversed labels, this is
    %needed bc the data retrieval using getCheckerBoard vs
    %getPrefixedEpochData is slightly different.
    
    %% Remove aftereffects? / remove bad muscles.
%     removeBadmuscles=0; %set this parameter outside.
    if removeBadmuscles==1
        if exist('musclesToRemove','var')
            normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData,musclesToRemove{1},musclesToRemove{2});
        else
            normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData);
        end
    end
end

Data={};%cell(7,1);
DataBySubj={};%cell(5,1); %intialize array to get data per epoch later, initialized size doesn't matter. the data array size will grow as needed.
% group=cell(5,1); %not used for now.
summFlag='nanmedian'; %how to summarize across subjects, used in plotCheckerboards to get correct data.


%% Pick muscles that you wanted to get the data from
wantedMuscles = {'TA','PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
newLabelPrefix = newLabelPrefix(contains(newLabelPrefix,wantedMuscles)); %if wantedMuscles is a full list, this will simply return the newLabelPrefix unchanged
n_muscles =length(wantedMuscles);
% muscleOrder={'TA','PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% newLabelPrefix2 = defineMuscleList(muscleOrder); %List of muscle
% newLabelPrefix2 = regexprep(newLabelPrefix2,'_s','s'); %removing the underscore "_"
% 
% for m=1:length(newLabelPrefix2)
%     
%     idx(m)=find(strcmp(newLabelPrefix,newLabelPrefix2(m))); %getting the index of the muscles
%     
% end
% wanted_Muscles= newLabelPrefix(sort(idx)); %It needs to be the muscle form the slow leg first
% newLabelPrefix= wanted_Muscles;
% n_muscles =length(muscleOrder); %number fo muscle in the experiment
%% Getting the data that we use for matrices
if contains(groupID,{'BATR','BATS'})
    ep=defineRegressorsDynamicsFeedback('nanmean'); % This is a hardcade (aka specific to this experiment) script to get the data from name epochs of interest
    epochOfInterest={'TM base','TM mid 1','PosShort_{early}',...
    'PosShort_{late}','Ramp','Optimal','Adaptation',...
    'Adaptation_{early}','TiedPostPos','OG2','NegShort_{early}','NegShort_{late}',...
    'Post1_{Early}','TMbase_{early}','Tied post Neg','Tied post ramp','OG base'}; % This line chooses the epochs that want to get data from
% epochOfInterest={'Post1_{Early}','TiedPostPos','TMmid2','Tied post ramp','Ramp','Tied post Neg'};
% epochOfInterest={'Ramp','Optimal'};
% epochOfInterest={'TM base','Ramp','Optimal','Tied post ramp'}; % This line chooses the epochs that want to get data from
    base=defineReferenceEpoch('TM base',ep); %pick the condition for baseline, used for remove bias.

elseif contains(groupID,{'AUF'}) %AUF study, load ep from Nirs_Automaticity repository
    if contains(normalizedGroupData.ID{1},'V03')
        intervention = true;
    else
        intervention = false;
    end
    ep = defineEpochNirs(intervention, 'nanmean'); %thie nan is used to summarize strides in an epoch
%     epochOfInterest = ep.Properties.ObsNames; %get all epochs 
    epochOfInterest = {'OGBase','TMBase','Adaptation_{SS}','NegShort_{la}','OGPost_{Early}','TMPost_{Early}','PosShort_{la}'};
    base = defineReferenceEpoch('TMBase',ep); %baseline condition when removing bias, here use TMBase to remove bias for Cs (TM conditions: adaptSS, negShortLa, PosShortLa)
end

% epochOfInterest = {'OGPost_{Early}'}; %{'NegShort_{la}','Adaptation_{SS}'};
% base = defineReferenceEpoch('OGBase',ep); %baseline condition when removing bias.

fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
ph=tight_subplot(1,length(epochOfInterest),[.03 .005],.04,.04);

flip=1; %1 for individual leg analysis and 2 for asymmetric (R-L), passed into groupAdaptationData.plotCheckerboards to get either individual or asymmetry data.
if flip==1
    n=2; %# of legs
    method='IndvLegs';
else
    n=1; %asym, only plot 1 set of muscles
    method='Asym';
end
fdr=.1; %currently not used?
C=[];

%set the removeBias outside.
% removeBias=1 %Flag for bias removal; to prep data for regression, should save the data without removeBias, bc will remove bias again later in model_fit_individual_muscles.m

for l=1:length(epochOfInterest)
    ep2=ep(strcmp(ep.Properties.ObsNames,epochOfInterest{l}),:); 
%     ep2=defineReferenceEpoch(epochOfInterest{l},ep); %locate the epoch of interests from the ep list. 
    %(this functions also renames the epoch of interest with prefix Ref,
    %which is not really needed in this case, only impact is adding a ref on the figure.)
    
    if removeBias==1
        [~,~,~,DataBySubj{l}]=normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),base,flip,summFlag); %plotting the data, this will return the data per subject: interval x muscles x epoch (1) x #subjects
        [~,~,~,Data{l}]=normalizedGroupData.getCheckerboardsData(newLabelPrefix,ep2,base,flip,summFlag); %get the summarzied data.

    else
        [~,~,~,DataBySubj{l}]=normalizedGroupData.plotCheckerboards(newLabelPrefix,ep2,fh,ph(1,l),[],flip,summFlag); %plotting the data
        [~,~,~,Data{l}]=normalizedGroupData.getCheckerboardsData(newLabelPrefix,ep2,[],flip,summFlag); %get the summarized data. The plotCheckerboards should have been able to return also the summarized data. This should be a feature improvements in labtool.
%         %getting the data from the plots; when there is only 1 subject,
%         this is the same as the line above.
    end
    
    C=[C reshape(Data{l}(:,end:-1:1),12*n_muscles*n,1)]; %Reshaping the data form the plot to a column vector per epoch of interest.
    %C size wil be 336 (12interval x 28 muscles) x #epochs of interest; if
    %doing asymmetry would be 168 (12interval x 14 muscles) x #epochs  
end

%%
%% Color definition
ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
% ex1=[1,0,0];
% ex2=[0,0,1];
cc=[0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350      0.0780    0.1840];

mid=ones(1,3);
N=100;
% gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];

%% Making the figure a bit prettier
fs=15; %font size
colormap(flipud(map)) %changing the color map to the one tha defined about, we flipup the matrix bc the code does L-R and we want R-L
% c=flipud('gray');
% colormap(flipud(gray))
set(gcf,'color','w'); %setting the background white
set(ph(:,1),'CLim',[-1 1]*1,'FontSize',fs); %making sure that the first plot color scheme goes [-1 1] and making the name of the labels larger
set(ph(:,2:end),'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
colorbar %Showing the colormap bar

%% Saving the data
if exist('saveResAndFigures','var') && saveResAndFigures
    % resDir = %[cd];% '/LTI models/']; %resDir is current directory
    if strcmpi(groupID,{'AUF'}) %group anaysis
        [resDir, ~, ~, ~, ~] = setupDataPath('ParamsForGroupPlot', 'V02', '', '');
        groupID = [groupID normalizedGroupData.ID{1}(end-2:end)];
    end
    save([resDir filesep groupID,'_',num2str(n_subjects), '_',method,'C',num2str(length(epochOfInterest)) ,'_ShortPertubations_RemovedBadMuscle_',num2str(removeBadmuscles), 'RemoveBias_',num2str(removeBias)], 'C','epochOfInterest')
    %% save graph
    saveas(fh, ['X:\Shuqi\NirsAutomaticityStudy\Data\GroupResults\Group22Sub\EMG\' groupID '_checker_allEpochs_groupMedian.png'])
end

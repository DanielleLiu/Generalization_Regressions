% This is a script to get the confidance interval of the step-by-step
% weights. The regressors are the same for each iteration and group 

%%
% % upload your path 
% main='/Users/dulcemariscal/Documents/GitHub/';
% % main='C:\Users\dum5\OneDrive - University of Pittsburgh\_BoxMigration\GitHub\';
% addpath(genpath([main, 'Generalization_Regressions']))
% addpath(genpath([main,'labTools']))
% addpath(genpath([main,'LongAdaptation']))
% addpath(genpath([main,'splitbelt-EMG-adaptation']))
% addpath(genpath([main,'EMG-LTI-SSM']))
% addpath(genpath([main,'matlab-linsys']))
% % rmpath(genpath([main,'PittSMLlab']))

%% Load data and Plot checkerboard for all conditions.
clear all; close all; clc;
saveFig = false;
groupID ={'AUF'};
EMGAllGroups=[];

if strcmp(groupID,'AUF')
    [groupDataPath, ~, ~, ~, visitnum] = setupDataPath('ParamsForGroupPlot', 'V04', '', '');
    rng(visitnum) %set seed for reproducibility for each session and to also have different samples per session.
    cd(groupDataPath)
    load('MusclesToRemove.mat')
end

for id=1:length(groupID)
    
    [normalizedGroupData, newLabelPrefix,n_subjects]=creatingGroupdataWnormalizedEMG(groupID{id});
    
    %% Removing bad muscles
    %This script make sure that we always remove the same muscle for the
    %different analysis; 
    normalizedGroupData= RemovingBadMuscleToSubj(normalizedGroupData, musclesToRemove{1},musclesToRemove{2});
    
    %% Getting the C values
    if contains(groupID,{'BAT'})
        % epochOfInterest={'TM base','TM mid 1','PosShort_{early}','PosShort_{late}','Ramp','Optimal','Adaptation','Adaptation_{early}','TiedPostPos','TMmid2','NegShort_{late}','Post1_{Early}','TMbase_{early}'};
        epochOfInterest={'TM base','NegShort_{late}','Ramp','Optimal'}; 
        %TMBase:-35 excluding last 5 (perhaps a mistake?); negshort_late: -10 excluding last 5
        %ramp: last 10 ignoring last 5; optimal: last40 excluding last 5; 
        ep=defineRegressorsDynamicsFeedback('nanmean');
    elseif contains(groupID, {'AUF'})
        if contains(normalizedGroupData.ID{1},'V03')
            intervention = true;
        else
            intervention = false;
        end
        ep = defineEpochNirs(intervention, 'nanmean');
        epochOfInterest = {'OGBase','TMBase','PosShort_{la}','Adaptation_{SS}','NegShort_{la}'};
    end
    
    flip=1;
    
    if flip==1
        n=2;
        method='IndvLegs';
    else
        n=1;
        method='Asym';
    end
    if id==1
        sub=1:n_subjects;
    else
        sub=n_subjects+1:13+n_subjects;
    end
    for s=1:n_subjects
        for iteIdx=1:length(epochOfInterest)
            ep2=defineReferenceEpoch(epochOfInterest{iteIdx},ep);
            adaptDataSubject = normalizedGroupData.adaptData{1, s};
            [~,~,~,Data{sub(s),iteIdx}]=adaptDataSubject.getCheckerboardsData(newLabelPrefix,ep2,[],flip);
            %the Data will be a #subj (14 if 1 group or 28 for 2 groups) x #epoch cell array, within each cell
            %is a 12x28 data matrix
        end
    end
    
    
    %%  Getting the step-by-step data
    % Adaptation epochs
    if contains(groupID,{'BAT'})
        strides=[-40 440 200];
        if contains(groupID{id},'TR') %for treadmill Post 1
            cond={'TM base','Adaptation','Post 1'}; %Conditions for this group
        else % for overground post 1
            cond={'OG base','Adaptation','Post 1'}; %Conditions for this group
        end
    elseif contains(groupID,{'AUF'})
        strides=[-40 890 440]; %Number per strides per condition
        cond={'OGBase','Adaptation','OGPost'}; %Conditions for this group
        %here use OGBase as baseline bc 1st post is OG
    end
    totalStrides = sum(abs(strides));
    post1Index = sum(abs(strides([1 2]))) + [1:5];
    
    exemptFirst=[1];
    exemptLast=[5]; %Strides needed
    names={};
    shortNames={};
    
    ep=defineEpochs(cond,cond,strides,exemptFirst,exemptLast,'nanmean',{'Base','Adapt','Post1'}); %epochs
    
    padWithNaNFlag=true; %If no enough steps fill with nan, let this on
    [dataEMG,labels,allDataEMG]=normalizedGroupData.getPrefixedEpochData(newLabelPrefix(end:-1:1),ep,padWithNaNFlag); %Getting the data in reverse order (sGLU first)
    
    %Flipping EMG (align left and right leg event, flip the 2nd leg's interval 7-12 first than 1-6 subintervals)
    for i=1:length(allDataEMG) 
        aux=reshape(allDataEMG{i},size(allDataEMG{i},1),size(labels,1),size(labels,2),size(allDataEMG{i},3));
        allDataEMG{i}=reshape(flipEMGdata(aux,2,3),size(aux,1),numel(labels),size(aux,4));
    end
    %allDataEMG is cell array of size #cond (e.g., Base, Adapt, Post1 for Dulce), for each cell: #strides per cond x 336 x #subjects
    %EMGdata is cell array concatenated together: 3D matrix: #totalStrides
    %x 336 x #subjects
    EMGdata=cell2mat(allDataEMG); %Getting EMG data per participants, 
    
    EMGAllGroups=cat(3, EMGAllGroups, EMGdata); %concatenate on the subject dimension
end

%% Bootstrapping
muscPhaseIdx=1:size(EMGdata,2); %
if any(contains(groupID,{'BAT'}))
    epochOfInterest={'TM base','NegShort_{late}','Ramp','Optimal'};
    context= find(strcmp(epochOfInterest,'Optimal')==1);
    % reactive=find(strcmp(epochOfInterest,'NegShort_{late}')==1);
    base=find(strcmp(epochOfInterest,'TM base')==1);
    reactive=find(strcmp(epochOfInterest,'Ramp')==1); %this is the old verion. new used negShort
    %Data TR, this includes 3 ep only, Base: last 40 baseline ignoring last 5.
    %first 440 of adaptation (ignore first), first 200 of post1 (ignoring 1st)
    TR=EMGAllGroups(:,:,1:12);
    TS=EMGAllGroups(:,:,13:24);
    n_toSampleC = 24; %24 total all used to get the C (regressors)
    n_toSampleY = 12; %12 per group
elseif contains(groupID,{'AUF'})
    base=find(strcmp(epochOfInterest,'TMBase')==1);
    %     reactive=find(strcmp(epochOfInterest,'PosShort_{la}')==1); %last 10 strides excluding the last 5 from Pos short ramp condition (already reached full slipt)
    context= find(strcmp(epochOfInterest,'Adaptation_{SS}')==1);%last 40 excluding last 5 in the split pos 18 condition (last split, 150 strides splits after 9 short exploration)
    reactive=find(strcmp(epochOfInterest,'NegShort_{la}')==1); %last 10 excluding last 5 in neg short (full neg split w/o ramp)
    %for AUF group, 1 group only same number to sample C and Y and same
    %data to reconstruct. For simplicity just use TR if only 1 group
    %exists.
    [n_toSampleC, n_toSampleY] = deal(n_subjects);
    TR = EMGAllGroups;
end

Cmuscles_unit=[];
bootstrap=1; %Do you want to run the loop (1 yes 0 No)
replacement=0; %do you want to do it with replacement of the data (y) (1 yes 0 No, if no will only loop through once.)
regre_Const=1; % To keep the regressors constants for both groups

if bootstrap
    if replacement
        n=2000; %number of iterations
    else
        n=1;
    end
    
    f = waitbar(0,'1','Name','Boostrapping Data',...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    
    setappdata(f,'canceling',0);
    
    unit=nan(12,2,28,n); %Creating nan matrices
    Cmuscles_inv=nan(12,2,28,n);
    Yhat=nan(12,totalStrides,28,n);
    temp4=nan(totalStrides,2,28,n);
    
    Xhat_TR = []; Yhat_TR=[]; dynamics_TR = []; YhatReactive_TR = []; YhatContext_TR = []; Ymuscles_TR_all = [];
    Ymuscles_TR_shifted = []; Yhat_TR_shifted = []; YhatReactive_TR_shifted = []; YhatContext_TR_shifted = [];
    R2 = struct();
    for iteIdx=1:n %loop for number of iterations
        
        DataBootInMatrix=[];
        groupDataByEpoch=[];
        DataBoot={};
        % Check for clicked Cancel button
        if getappdata(f,'canceling')
            break
        end
        % Update waitbar and message
        ww=waitbar(iteIdx/n,f,['Iteration ' num2str(iteIdx)]);
        
        if  regre_Const %out of 24 bc both group combined to get the C
            IdxToSampleC=datasample(1:n_toSampleC,n_toSampleC,'Replace',false);
        else
            IdxToSampleC=datasample(1:n_toSampleC,n_toSampleC,'Replace',true);
        end
        
        if replacement %doing the bootstrap with replacement
            %sample out of 12 bc genereating index per group, Y is used per
            %group.
            IdxToSampleY=datasample(1:n_toSampleY,n_toSampleY,'Replace',true);
        else
            IdxToSampleY=datasample(1:n_toSampleY,n_toSampleY,'Replace',false);
        end
        
        DataBoot=Data(IdxToSampleC,:); %Subject pick at each loop for the Cs;
        %the Data will be a #subj (14 if 1 group or 28 for 2 groups) x #epoch cell array, within each cell
        %is a 12x28 data matrix
        
        %This loop is to compute the constrant on our regressions
        if iteIdx==1 %get it only once, same for all iterations.
            for c=1:length(epochOfInterest)
                for s=1:length(IdxToSampleC) %transform DataBoot from cellarray that includes #subjects of 2D matrix to a 3D matrix
                    DataBootInMatrix(:,:,s)=DataBoot{s,c}; %temp size would be 12x28x24 (s=24), last dimension =#subjects
                end
                groupDataByEpoch{1,c}=nanmedian(DataBootInMatrix,3); %median over subjects
                tt=groupDataByEpoch{1,c}(:,end:-1:1); %get the data with muscle order (2nd dimension reversed)
                groupDataByEpoch{1,c}=tt;
                groupDataByEpoch{1,c}=reshape(groupDataByEpoch{1,c},14*2*12,1); %reshping the data for the C values, reshape to column vectors for each epoch of interest
                %x final size: cell array of size epoch of interest, each
                %cell contains a 336x1 column vector which is the median of
                %subjects and each subject's data is mean of the epoch.
            end
            
            
            groupDataByEpoch=cell2mat(groupDataByEpoch)';
            groupDataByEpoch=groupDataByEpoch';
            
            %C values that we are using for the regressions,
            %Always remove bias, this base is defined by 
            %ep = defineRegressorsDynamicsFeedback or defineEpochNirs
            %from the beginning of the code. For BAT: seems to be last 35
            %excluding last 5 (FIXME: perhaps a mistake); For AUF: last 40
            %excluding last 5. usually average of the epoch.
            Cmuscles=[groupDataByEpoch(:,reactive) groupDataByEpoch(:,context)]- groupDataByEpoch(:,base);
           
            %reorganize the data to be separatend by muscle
            Cmuscles=reshape(Cmuscles',2,12,28); 
            
        end
        %Picking the data muscles that we want and participants
        Y_TR=TR(:,muscPhaseIdx,IdxToSampleY); %groupIdx could be sampled with or without replacement
        
        %removing the bias for group
        Y_TR=nanmedian(Y_TR,3); %getting the median of the group 
        %FIXME: this is removing different bias for Y and for X. not sure
        %why. Oringally 5:30 average over 5-30 strides from the last 40 (excluding last 5); X
        %bias is average of the last 35 (excluding last 5) by definition of
        %defineRegressorsDynamicsFeedback. Fixed to be 6-40 for BAT data
        %for AUF: should always do avg of last 40 excluding last 5.
        %FIXME: notice also we are removing OGbaseline for all periods in
        %TS group and the regressor used TMBaseline removal
        bias=nanmean(Y_TR(1:40,:,:)) ; %estimating the gorup baseline, assuming first 40 strides are TMBase.
        Y_TR=Y_TR-bias; %removing the bias from the data
           
        %reorganize the data to be separatend by muscle
        Ymuscles_TR=reshape(Y_TR,[],12,28);
        
        if any(contains(groupID,{'BATS'}))
            Y_TS=TS(:,muscPhaseIdx,IdxToSampleY);
            Y_TS=nanmedian(Y_TS,3); %getting the median of the group 
            %FIXME: here similarly stride index seems weird.
            bias=nanmean(Y_TS(6:40,:,:)) ; %estimating the group baseline
            Y_TS=Y_TS-bias; %removing the bias from the data
            Ymuscles_TS=reshape(Y_TS,[],12,28);
        end     
        
        %Linear regression individual muscles
        reconstruction_indv=[];
        data=[];
        C_indv=[];
        X_indv=[];
        reconstruction_indv_shifted=[];yhat_context_shifted=[]; yhat_reactive_shifted=[]; data_shifted=[];
        yhat_reactive=[];yhat_context=[];
        for i=1:size(Ymuscles_TR,3) %loop for individual muscle fit
            
           % Getting the inverse of the regressors
%            if l==1
            unit=Cmuscles(:,:,i)'./vecnorm(Cmuscles(:,:,i)');
            if iteIdx==1
                Cmuscles_unit=[Cmuscles_unit;unit];
            end
            Cmuscles_inv=pinv(unit'); %geeting the inverse of the constant


             %%% TR
            Xhat_TR(:,:,i,iteIdx) =Cmuscles_inv'*Ymuscles_TR(:,:,i)'; %x= y/C
            Yhat_TR(:,:,i,iteIdx) =  unit* Xhat_TR(:,:,i,iteIdx) ; %Estimated Y with the constants, 12 x stride x 28 x iterations
            dynamics_TR(:,:,i,iteIdx)=Xhat_TR(:,:,i,iteIdx)'; %step-by-step dynamics, e.g., weights. 2 x #strides x 28 x iterations
            %unit is reactive first then context.
            YhatReactiveCurr = unit(:,1)* Xhat_TR(1,:,i,iteIdx) ; %Estimated Y reactive, 1 (reactive weight) x stride x 28 (muscles) x iterations
            YhatContextCurr =  unit(:,2)* Xhat_TR(2,:,i,iteIdx) ; %Estimated Y context, 1 (context weight) x stride x 28 x iterations
            yhat_context=[yhat_context; YhatContextCurr];
            yhat_reactive=[yhat_reactive; YhatReactiveCurr];
            
            r = abs(nanmin(YhatReactiveCurr)); %reactive min
            c = abs(nanmin(YhatContextCurr)); %contextual min
            
            reconstruction_indv_shifted =[reconstruction_indv_shifted ;  Yhat_TR(:,:,i,iteIdx)+r+c]; % Concatenating the data reconstructed
            yhat_context_shifted=[yhat_context_shifted; YhatContextCurr+c];
            yhat_reactive_shifted=[yhat_reactive_shifted; YhatReactiveCurr+r];
            data_shifted =[ data_shifted ; Ymuscles_TR(:,:,i)'+c+r];  % Concatenating the data
%             Ymuscles_TR_shifted(:,:,i,iteIdx) = Ymuscles_TR(:,:,i)' + c +r; %12 x stride x 28
%             Yhat_TR_shifted(:,:,i,iteIdx)=  Yhat_TR(:,:,i,iteIdx) + c +r ; %Estimated Y with the constants, 12 x stride x 28 x iterations
%             YhatReactive_TR_shifted(:,:,i,iteIdx) = YhatReactive_TR(:,:,i,iteIdx) + r; %Estimated Y reactive, 1 (reactive weight) x stride x 28 (muscles) x iterations
%             YhatContext_TR_shifted(:,:,i,iteIdx) =  YhatContext_TR(:,:,i,iteIdx) + c; %Estimated Y context, 1 (context weight) x stride x 28 x iterations
            
            if any(contains(groupID,{'BATS'}))
                %%% TS
                Xhat_TS(:,:,i,iteIdx) =Cmuscles_inv'*Ymuscles_TS(:,:,i)'; %x= y/C
                Yhat_TS(:,:,i,iteIdx)=  unit* Xhat_TS(:,:,i,iteIdx) ; %Estimated Y with the constants
                dynamics_TS(:,:,i,iteIdx)=Xhat_TS(:,:,i,iteIdx)'; %step-by-step dynamics
            end 
        end
        data = reshape(nanmean(Ymuscles_TR(post1Index,:,:),1),[],1);
        yhat = reshape(nanmean(Yhat_TR(:,post1Index,:,iteIdx),2),[],1);
        yhat_context = nanmean(yhat_context(:,post1Index),2);
        yhat_reactive = nanmean(yhat_reactive(:,post1Index),2);
        
        R2.shifted0.relativeToMean0(1,iteIdx) = my_Rsquared_coeff(data,yhat,false);
        R2.shifted0.relativeToMean0(2,iteIdx) = my_Rsquared_coeff(data,yhat_context,false);
        R2.shifted0.relativeToMean0(3,iteIdx) = my_Rsquared_coeff(data,yhat_reactive,false);

        R2.shifted0.relativeToMean1(1,iteIdx) = my_Rsquared_coeff(data,yhat,true);
        R2.shifted0.relativeToMean1(2,iteIdx) = my_Rsquared_coeff(data,yhat_context,true);
        R2.shifted0.relativeToMean1(3,iteIdx) = my_Rsquared_coeff(data,yhat_reactive,true);

        data_shifted = nanmean(data_shifted(:,post1Index),2);
        reconstruction_indv_shifted = nanmean(reconstruction_indv_shifted(:,post1Index),2);
        yhat_context_shifted = nanmean(yhat_context_shifted(:,post1Index),2);
        yhat_reactive_shifted = nanmean(yhat_reactive_shifted(:,post1Index),2);
        
        R2.shifted1.relativeToMean0(1,iteIdx) = my_Rsquared_coeff(data_shifted,reconstruction_indv_shifted,false);
        R2.shifted1.relativeToMean0(2,iteIdx) = my_Rsquared_coeff(data_shifted,yhat_context_shifted,false);
        R2.shifted1.relativeToMean0(3,iteIdx) = my_Rsquared_coeff(data_shifted,yhat_reactive_shifted,false);

        R2.shifted1.relativeToMean1(1,iteIdx) = my_Rsquared_coeff(data_shifted,reconstruction_indv_shifted,true);
        R2.shifted1.relativeToMean1(2,iteIdx) = my_Rsquared_coeff(data_shifted,yhat_context_shifted,true);
        R2.shifted1.relativeToMean1(3,iteIdx) = my_Rsquared_coeff(data_shifted,yhat_reactive_shifted,true);
         
    end
    close(ww)
    
end

delete(f)

% %in case previosu error couldn't close wait bar
% F = findall(0,'type','figure','tag','TMWWaitbar')
% delete(F)
R2.labels = {'total','contextual','reactive'};

%% compute R2 of the fit
iteIdx=n;
for strideToPlot = 1:1370
    yactualCurr = reshape(squeeze(Ymuscles_TR(strideToPlot,:,:)),[],1);%FIXME: this is incorrect, only have the last iteration
    yhatCurr = reshape(squeeze(Yhat_TR(:,strideToPlot,:,iteIdx)),[],1); 
    R2Example(strideToPlot) = my_Rsquared_coeff(yactualCurr,yhatCurr, true);
end
figure(); plot(R2Example); title('R^2 of the last iteration');

%% Plot regressors checkerboards
strideToPlot = 931:935; %[40 890 440], only make sense to check post-adapt bc baseline removal removed OGBase from the whole exp which is only appropriate for OGPost.
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle 
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);

figure
subplot(1,5,1)
imagesc((reshape(Cmuscles_unit(:,1),12,28)'))
title('Reactive')
fs = 10;
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
set(gca,'CLim',[-1 1]*1,'FontSize',fs); %making sure that the first plot color scheme goes [-1 1] and making the name of the labels larger

subplot(1,5,2)
imagesc((reshape(Cmuscles_unit(:,2),12,28)'))
title('Contextual')

subplot(1,5,3)
imagesc(squeeze(nanmean(Yhat_TR(:,strideToPlot,:,iteIdx),2))') %avg over strides
title('Yhat')

subplot(1,5,4)
imagesc(squeeze(nanmean(Ymuscles_TR(strideToPlot,:,:),1))') %avg overstrides
title('Y')

y_original = reshape(Ymuscles_TR,size(Ymuscles_TR,1),[]) + bias;
y_original = reshape(y_original,size(y_original,1),12,28);
subplot(1,5,5)
imagesc(squeeze(nanmean(y_original(strideToPlot,:,:),1))') %avg overstrides
title('YOriginal')

% %get Y before removing bias.
% bias=nanmean(Y_TR(1:40,:,:)) ; %estimating the gorup baseline, assuming first 40 strides are TMBase.
% Y_TR=Y_TR-bias; %removing the bias from the data
% 
% %reorganize the data to be separatend by muscle
% Ymuscles_TR=reshape(Y_TR,[],12,28);
        
color4checkerboards


% Making the figure a bit prettier
fs=14; %font size
colormap(flipud(map)) %changing the color map to the one tha defined about, we flipup the matrix bc the code does L-R and we want R-L
% c=flipud('gray');
% colormap(flipud(gray))
% set(gcf,'color','w'); %setting the background white
set(gca,'YTickLabels',{},'CLim',[-1 1]*1,'FontSize',fs);
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
colorbar %Showing the colormap bar

% compute R2 of the fit

sgtitle(['Last Iteration Fit R^2 = ' num2str(nanmean(R2Example(strideToPlot))) ' Stride=' num2str(strideToPlot)])


%% plot the CI of the bootstrapped results.
postBoostrapAvg = []; %reactive 1st row, context next row
titles = {'Reactive','Contextual'};
for cIdx=1:2 %reactive, then context, 
    f = figure('units','normalized','outerposition',[0 0 1 0.5]);%('Position', get(0, 'Screensize'));
    for mIdx = 1:size(dynamics_TR,3) %for each muscle  muscle orders start from sGLU
       postData = nanmean(squeeze(dynamics_TR(post1Index,cIdx,mIdx,:)));
       PlotHelper.plotCI(mIdx, postData,'b',ytl{mIdx},false, true) %
       postBoostrapAvg(cIdx,mIdx) = nanmean(postData);
    end
    title([titles{cIdx} ' (Mean, 95% from bootstrap)'])
    xlim([0 size(dynamics_TR,3)+1])
    xticks(1:size(dynamics_TR,3))
    xticklabels(ytl)
    xtickangle(45)
    PlotHelper.tightMargin(gca)
    if saveFig && any(contains(groupID,{'AUF'}))
        resDir = ['X:\Shuqi\NirsAutomaticityStudy\Data\GroupResults\Group22Sub\EMG\V02' filesep]
        exportgraphics(f,[resDir 'Btstrp' titles{cIdx} '.png'],'Resolution',300);
    end
end

postBootstrapData.avg = postBoostrapAvg;
postBootstrapData.muscleLabels = ytl;
postBootstrapData.stateLabels = titles;
% save('GroupBootstrapDataOGPost.mat','postBootstrapData')

%% plot VAF from bootstrap
saveFigPath = 'X:\Shuqi\NirsAutomaticityStudy\Data\GroupResults\Group22Sub\EMG\V04\Bootstrp2000_';
subIds = sprintfc('%d',1:n);
for shifted = 0:1
    for relativeToMean = 0:1
        eval(['dataToPlot = R2.shifted' num2str(shifted) '.relativeToMean' num2str(relativeToMean) ';'])
        PlotHelper.barPlotWithIndiv(dataToPlot,subIds,{'Total','Context','Reactive'},'VAF',['OGPost VAF Shifted', num2str(shifted) ' RelativeToMean' num2str(relativeToMean)], ...
            saveFig,[saveFigPath 'VAF_Shifted_', num2str(shifted) '_RelativeToMean_' num2str(relativeToMean) '_boostrap'],[],[],[],false);
        legend('NumColumns',2)
    end
end

%% save the boot strap data.
if any(contains(groupID,{'BAT'}))
    save([groupID{1}(1:3),'_',num2str(sub(end)),'_iteration_', num2str(iteIdx),'_Individual_muscles_C_constant_Adaptation_per_group'],'dynamics_TR','dynamics_TS','-v7.3')
% save([groupID{1}(1:3),'_',num2str(24),'_iteration_', num2str(l),'_Individual_muscles_C_constant_Adaptation'],'dynamics_TR','-v7.3')
% save([groupID,'_',num2str(n_subjects),'_iteration_', num2str(n),'_Individual_muscles'],'dynamics','Yhat','Ymuscles','groupID','-v7.3')
elseif any(contains(groupID,{'AUF'}))
    save([groupID{1}(1:3),normalizedGroupData.ID{1}(end-2:end),'_',num2str(sub(end)),'_iteration_', num2str(iteIdx),'_Individual_muscles_C_constant_PostAdaptation_per_group'],'dynamics_TR','R2','-v7.3')
end

return %stop the script here, below are different plotting options.

%%
load('musclesLabels.mat')
% load BAT_24_iteration_2000_Individual_muscles.mat
% load('BATS_12_iteration_2000_Individual_muscles.mat')
OG=dynamics_TS;
% % load('BATR_12_iteration_2000_Individual_muscles.mat')
TM=dynamics_TR;


% load('musclesLabels.mat')
% load('BATR_indv_muscles.mat')
% TM_2=X2asym;
% load('BATS_indv_muscles.mat')
% OG_2=X2asym;
%%
load('NCM2023_Treadmill.mat')
TM_2=X2asym;
TM_2(1,:,:)=-TM_2(1,:,:);
load('NCM2023_OG.mat')
OG_2=X2asym;
OG_2(1,:,:)=-OG_2(1,:,:);
%% Group data plotting 
clrMap = colorcube(28*3);
muscles=[1:14 1:14];
g=[14:-1:1 14:-1:1];
ff=[1:14 1:14;15:28 15:28];
range=481:485;


for i=1:28
DataBootInMatrix{i,1}=labels(i).Data(1:end-1);
end

%contextual 
% label2={'fLG','fMG'}; %Yes TM and Yes OG 
% label2={'sTFL','sHIP','sVM','sPER','sSOL','sSEMT','sSEMB','sLG','sMG'}; %Yes TM; No OG
% label2={'fGLU','fTFL','fHIP','fRF','sRF','fVL','sVL','fVM','fSEMT','fSEMB','sTA','sBF','fBF'}; %No TM and NO OG 
% label2={'sGLU','fPER','fSOL','fTA'}; %no TM and Yes OG 


%reactive 
% label2={'fGLU','sTA','fTA','sPER','fPER','sSOL','fSOL','sLG','fLG','fMG','sMG','fBF','sBF',...
% 	'sSEMB','sSEMT','fVM','fVL','fRF','fHIP','sHIP','fTFL'}; %Yes TM and Yes OG 
% label2={'sGLU','fSEMB','fSEMT'};  %Yes TM; No OG
% label2={'sTFL'};  %No TM and NO OG 
label2={'sVM','sVL','sRF'}; %NO TM and YES OG 
t=[];
t(1,:)=find(contains(DataBootInMatrix,label2));
%%
for dyn=1%
    figure()
    hold on
    DataBootInMatrix=[];
    groupDataByEpoch=[];
    
    for m=t
%         figure(ff(dyn,m))
%         hold on
        groupDataByEpoch=squeeze(TM(:,dyn,m,:));
        y=squeeze(OG(:,dyn,m,:));
        
        x_mean=nanmean(groupDataByEpoch(range,:),'all');
        y_mean=nanmean(y(range,:),'all');
        
        centers=[x_mean  y_mean];
        
        P_x = prctile(nanmean(groupDataByEpoch(range,:),1)',[2.5 97.5],"all");
        P_y = prctile(nanmean(y(range,:),1)',[2.5 97.5],"all");
        
        llc=[P_x(1), P_y(1)];
        
        CIrng(1)=P_x(2)-P_x(1);
        CIrng(2)=P_y(2)-P_y(1);
        
        x0=x_mean; % x0,y0 ellipse centre coordinates
        y0=y_mean;
 
        
        text(x0+.02,y0,{labels(m).Data(1:end-1)})
        
        if P_y(1)<0 &&  P_y(2)>0 && P_x(1)<0 &&  P_x(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',"-",'LineWidth',1);
        elseif P_y(1)<0 &&  P_y(2)>0 
             rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',"-",'LineWidth',1);
        elseif P_x(1)<0 &&  P_x(2)>0
               rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','-','LineWidth',1);
        else 
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineWidth',1);
        end

        
        plot(centers(1), centers(2), 'o', 'MarkerFaceColor', clrMap(m+3,:), 'MarkerSize',10, 'LineWidth', 1,'MarkerEdgeColor',clrMap(m+3,:))%clrMap(m+3,:))
        xlabel('Treadmill')
        ylabel('Overground')

    end
    xlabel('Treadmill')
    ylabel('Overground')
%     
    if dyn==1
        title('Reactive')
        xx=-.5:0.1:2.1;
        plot(xx,xx,'r')
%         xlim([-.5 2.5])
%         ylim([-.5 2.5])
        yline(0)
        xline(0)
    else
        title('Contextual')
        yline(0)
        xline(0)
        xlim([-.6 1.6])
        ylim([-.6 1.6])
    end
    
    set(gcf,'color','w')
  
end

%% Individual muscle data plotting
clrMap = colorcube(28*3);
muscles=[1:14 1:14];
g=[14:-1:1 14:-1:1];
ff=[1:14 1:14;15:28 15:28];
range=481:485; %Post-adapt 481:485 Early-adapt 41:45 Late adapt 440:480
% range=
for dyn=1
    
    
    DataBootInMatrix=[];
    groupDataByEpoch=[];
    
    for m=1:28
        figure(ff(dyn,m))
        hold on
        groupDataByEpoch=squeeze(TM(:,dyn,m,:));
        y=squeeze(OG(:,dyn,m,:));
        
        x_mean=nanmean(groupDataByEpoch(range,:),'all');
        y_mean=nanmean(y(range,:),'all');
        
        centers=[x_mean  y_mean];
        
        P_x = prctile(nanmean(groupDataByEpoch(range,:),1)',[2.5 97.5],"all");
        P_y = prctile(nanmean(y(range,:),1)',[2.5 97.5],"all");
        
        llc=[P_x(1), P_y(1)];
        
        CIrng(1)=P_x(2)-P_x(1);
        CIrng(2)=P_y(2)-P_y(1);
        
        x0=x_mean; % x0,y0 ellipse centre coordinates
        y0=y_mean;
        
        
        text(x0+.02,y0,{labels(m).Data(1:end-1)})
        
        %         if P_y(1)<0 &&  P_y(2)>0 || P_x(1)<0 &&  P_x(2)>0
        if P_y(1)<0 &&  P_y(2)>0 && P_x(1)<0 &&  P_x(2)>0
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',":",'LineWidth',1);
        elseif P_y(1)<0 &&  P_y(2)>0 
             rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle',"-.",'LineWidth',1);
        elseif P_x(1)<0 &&  P_x(2)>0
               rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineStyle','--','LineWidth',1);
        else 
            rectangle('Position',[llc,CIrng],'Curvature',[1,1],'EdgeColor',clrMap(m+3,:),'LineWidth',2);
        end
        
        
        plot(centers(1), centers(2), 'o', 'MarkerFaceColor', clrMap(m+3,:), 'MarkerSize',10, 'LineWidth', 1,'MarkerEdgeColor',clrMap(m+3,:))%clrMap(m+3,:))
        xlabel('Treadmill')
        ylabel('Overground')
        %
%         if m<15
%             Li{1}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'b');
%             text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
%         else
%             Li{2}=scatter(nanmean(TM_2(dyn,range,m)),nanmean(OG_2(dyn,range,m)),100,"filled",'MarkerFaceColor', 'r')  ;
%             text(nanmean(TM_2(dyn,range,m))+.02,nanmean(OG_2(dyn,range,m)),{labels(m).Data(1:end-1)})
%         end
        
        if dyn==1
            title('Reactive')
            xx=-.5:0.1:2.1;
            plot(xx,xx,'r')
            
%             xlim([-.5 2])
%             ylim([-.5 2])
            yline(0)
            xline(0)
        else
            title('Contextual')
            yline(0)
            xline(0)
%             xlim([-1 1])
%             ylim([-1 1])
        end
        xlabel('Treadmill')
        ylabel('Overground')
        set(gcf,'color','w')
    end
    
end
%% 
%%Time Courses
colors=[0 0.4470 0.7410;0.8500 0.3250 0.0980];
adaptation=0
if adaptation==1
    range=1:480;
else
    range=481:680;
end

for m=1:28
    
    DataBootInMatrix=[];
    groupDataByEpoch=[];
    figure
    for dyn=1:2
        groupDataByEpoch=squeeze(TM(:,dyn,m,:));
        y=squeeze(OG(:,dyn,m,:));
        
        x_mean=nanmean(groupDataByEpoch(range,:),2);
        y_mean=nanmean(y(range,:),2);
        
        P_x = prctile(groupDataByEpoch(range,:)',[2.5 97.5],1);
        P_y = prctile(y(range,:)',[2.5 97.5],1);
        %     x = 1:numel(y);
        %         x=index{i};
        %         std_dev = nanstd(d{i},1);
        
        subplot(2,1,1)
        hold on
        curve1 = P_y(1,:)';
        curve2 = P_y(2,:)';
        y2 = [1:size(x_mean,1), fliplr(1:size(x_mean,1))];
        inBetween = [curve1', fliplr(curve2')];
        fill(y2, inBetween,colors(dyn,:),'FaceAlpha',0.3,'EdgeColor','none')
        hold on;
        ylabel('W')
        xlabel('Strides')
        Li{dyn}=plot(1:size(x_mean,1),y_mean,'LineWidth', 2,'Color',colors(dyn,:));
        title(['Overground' ,{labels(m).Data(1:end-1)}])
        yline(0)
        
        subplot(2,1,2)
        hold on
        curve1 = P_x(1,:)';
        curve2 = P_x(2,:)';
        y2 = [1:size(x_mean,1), fliplr(1:size(x_mean,1))];
        inBetween = [curve1', fliplr(curve2')];
        fill(y2, inBetween,colors(dyn,:),'FaceAlpha',0.3,'EdgeColor','none')
        hold on;
        ylabel('W')
        xlabel('Strides')
        title(['Treadmill' ,{labels(m).Data(1:end-1)}])
        Li{dyn}=plot(1:size(x_mean,1),x_mean, 'LineWidth', 2,'Color',colors(dyn,:));
        yline(0)
        axis tight
    end
    
    legend([Li{:}],[{'Reactive';'Contextual'}])
    set(gcf,'color','w')
end
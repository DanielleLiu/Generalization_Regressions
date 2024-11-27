%Model fit per muscle

%% Adding path 

% % upload your path
% main='/Users/dulcemariscal/Documents/GitHub/';
% % main='C:\Users\dum5\OneDrive - University of Pittsburgh\_BoxMigration\GitHub\';
% addpath(genpath([main, 'Generalization_Regressions']))
% addpath(genpath([main,'labTools']))
% % addpath(genpath([main,'LongAdaptation']))
% addpath(genpath([main,'splitbelt-EMG-adaptation']))
% addpath(genpath([main,'EMG-LTI-SSM']))
% addpath(genpath([main,'matlab-linsys']))
% % rmpath(genpath([main,'PittSMLlab']))

%% Load real data:
clear all;clc;close all
% Free model - Linear regression - Asymmetry with baseline
%% This is just the saved data - Update accrodingly 
groupID='BAT';
adapt = true; %false for OGPost.
groupAnalysis = false;
cd('/Volumes/Research/Shuqi/Grant/2023Summer/Data/');

if contains(groupID,'BAT')
% groupID='C3'
%EMG Data 
% fname = ['dynamicsData_', groupID, '_subj_3_RemoveBadMuscles_0.h5']
if adapt
    fname='dynamicsData_BATS_subj_12_RemoveBadMuscles_1_Adaptation.h5'
else
    fname ='dynamicsData_BATS_subj_12_RemoveBadMuscles_1_PostAdaptation.h5';
end
%%%%%%%% Regressors
% load([ groupID,'_1_IndvLegsC4_RemovedBadMuscle_0RemoveBias_1.mat'])
%%%%%%%% Regressors
%BATR
% load BATR_12_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BATR_12_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

%BATS
% load BATS_12_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemovBias_0.mat

% ALL 24 participants
% load BAT_24_AsymC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat
% load BAT_24_IndvLegsC17_ShortPertubations_RemovedBadMuscle_1RemoveBias_0_PosteriorMuscles.mat
elseif contains(groupID,'AUF')
    if contains(groupID,'V04')
        fname = '21';
        [groupDataPath, ~, ~, ~, ~] = setupDataPath('ParamsForGroupPlot', 'V04', '', '');
%         cd(groupDataPath)
    else
        fname = '22';
        [groupDataPath, ~, ~, ~, ~] = setupDataPath('ParamsForGroupPlot', 'V02', '', '');
%         cd(groupDataPath)
    end
    if adapt
        fname = ['dynamicsData_AUF_subj_' fname '_RemoveBadMuscles_1_Adapt_TMBase.h5'];
    else
        fname = ['dynamicsData_AUF_subj_' fname '_RemoveBadMuscles_1_OGPost_OGBase.h5']; %h5 file name
    end
end
epochIdxToPlot = 41:45; %for earlyPost or earlyAda
%%
muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);
fs=14;
%% Getting the data
EMGdata=h5read(fname,'/EMGdata');
binwith=10;
[Y,~,U,Ubreaks,Ysum,Yindv,labels,C,Cinv]=groupDataToMatrixForm_Update(1:size(EMGdata,3),fname,0);

% color=colormap(jet(size(Yindv,3)));
poster_colors;
colorOrder=[p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1];p_red; p_orange; p_fade_green; p_fade_blue; p_plum; p_green; p_blue; p_fade_red; p_lime; p_yellow; [0 0 0];[0 1 1]];
colorOrder=[ colorOrder; colorOrder;colorOrder];

%% Organizing the data

subID=h5read(fname,'/SubID')';
subID=deblank(subID); %Remove random null values added during the saving of the matrix

Uf=[U;ones(size(U))];
removebaseline=1; %Remove bias flag
c_constant=0;
showPlot=1;
figure

data_perSubject=[]; yhatContext_perSubject=[]; C_perSubject=[]; yhat_perSubject=[]; %336 x strides x subject.
yhatReactive_perSubject=[]; 
yhatContext_perSubject_shifted=[]; yhatReactive_perSubject_shifted = []; yhat_perSubject_shifted=[]; data_perSubject_shifted=[];

if groupAnalysis
    n_sub = 1;
else
    n_sub = size(Yindv,3);
end
mdl_perSubject=cell(n_sub,28);

for subj=1:n_sub
    if groupAnalysis
        Yasym = Y;
    else
        Yasym=Yindv(:,:,subj);
    end
    Ymodel=Yasym';  % Transpose the EMG data to match equations

    if c_constant==1
        %Pick the conditions that we are going to use for the regression
        %SL: 08/16/2023: the following is obsolete, should load h5 file
        %only right now.
        if contains(groupID,'BAT')
            load BAT_24_IndvLegsC16_ShortPertubations_RemovedBadMuscle_1.mat
            reactive=find(strcmp(epochOfInterest,'Ramp')==1);
        elseif contains(groupID,'AUFV02')
            load AUFV02_22_IndvLegsC18_ShortPertubations_RemovedBadMuscle_0RemoveBias_0.mat
        end
    else 
        if ~groupAnalysis
            C=Cinv(:,:,subj)';
        end
        %Pick the conditions that we are going to use for the regression
        epochOfInterest=h5read(fname,'/Epochs')';
        epochOfInterest=deblank(epochOfInterest); %Remove random null values added during the saving of the matrix

    end
    if contains(groupID,'BAT')
        context= find(strcmp(strtrim(epochOfInterest),'Optimal')==1);
        reactive2=find(strcmp(epochOfInterest,'NegShort_{early}')==1);
        reactive=find(strcmp(epochOfInterest,'Ramp')==1);
        tmbase=find(strcmp(epochOfInterest,'TM base')==1);
    elseif contains(groupID,'AUF')
        reactive=find(strcmp(epochOfInterest,'PosShort_{la}')==1); %last 10 strides excluding the last 5 from Pos short ramp condition (already reached full slipt)
        context= find(strcmp(epochOfInterest,'Adaptation_{SS}')==1);%last 40 excluding last 5 in the split pos 18 condition (last split, 150 strides splits after 9 short exploration)
        reactive2=find(strcmp(epochOfInterest,'NegShort_{la}')==1); %last 10 excluding last 5 in neg short (full neg split w/o ramp)
        tmbase=find(strcmp(epochOfInterest,'TMBase')==1);
    end
    if adapt %adaptation use pos short as reactive
        Casym=[C(:,reactive) C(:,context)]; % EMGreactive and EMGcontext
    else
        Casym=[C(:,reactive2) C(:,context)]; % EMGreactive and EMGcontext
    end




if removebaseline==1
    bias=nanmean(Yasym(1:40,:)); %Computing the bias
    
    Casym=Casym-C(:,tmbase); %Removing bias from the regressors
%      Casym=Casym-bias'; %Removing bias from the regressors
    Ymodel=Ymodel-bias'; %removing bias from the data
end
    %Organizing the data per muscle
    Cmuscles=reshape(Casym',2,12,size(Yasym,2)/12); %muscles for regressors 2 number of regressors, 12 number of phase of the gait cycle
    Ymuscles=reshape(Ymodel(:,1:size(Ymodel,2))',size(Ymodel,2),12,size(Yasym,2)/12); %data

% Linear regression individual muscles
reconstruction_indv=[];
data=[];
C_indv=[];
yhatContext = [];
yhatReactive = [];
reconstruction_indv_shifted=[];
yhat_reactive_shifted=[];
yhat_context_shifted=[];
data_shifted = [];
% X_indv=[];

for i=1:28
    unit(:,:,i)=Cmuscles(:,:,i)'./vecnorm(Cmuscles(:,:,i)'); %Getting the unit vector of the regressors
    if all(isnan(unit(:,:,i)))
        temp(:,:,i)=(unit(:,:,i)); %if nan can't get inverse, carry the nan alone.
    else 
        temp(:,:,i)=pinv(unit(:,:,i)'); % Getting the inverse, if data is not nan.
    end
    X2asym(:,:,i) =temp(:,:,i)'*Ymuscles(:,:,i)'; %x= y/C
    Y2asym(:,:,i)=  unit(:,:,i)* X2asym(:,:,i) ; %yhat
    temp2(:,:,i)=X2asym(:,:,i)'; % transposing the dynamics vector
    model{i}.C=unit(:,:,i); %saving the regressors. Yhis is necesary for the plotting funciton
    
    reconstruction_indv =[ reconstruction_indv ; Y2asym(:,:,i)]; % Concatenating the data reconstructed
    data =[ data ; Ymuscles(:,:,i)'];  % Concatenating the data
    C_indv=[C_indv;unit(:,:,i)];  % Concatenating the regressors for each muscles
    yContextCurr = unit(:,2,i)* X2asym(2,:,i);
    yReactiveCurr = unit(:,1,i)* X2asym(1,:,i);
    yhatContext = [yhatContext; yContextCurr ]; %the contextual component of the yhat
    yhatReactive = [yhatReactive; yReactiveCurr]; %the contextual component of the yhat
%     X_indv=temp2(:,:,i);  % Concatenating the dynamics

    %shift to be all positive
%     r=abs(nanmin(unit(:,2,i)));
%     c=abs(nanmin(unit(:,1,i)));
    r = abs(nanmin(yReactiveCurr)); %reactive min
    c = abs(nanmin(yContextCurr)); %contextual min
    
    reconstruction_indv_shifted =[ reconstruction_indv_shifted ; Y2asym(:,:,i)+c+r]; % Concatenating the data reconstructed
    yhat_context_shifted=[yhat_context_shifted; yContextCurr+c];
    yhat_reactive_shifted=[yhat_reactive_shifted; yReactiveCurr+r];
    data_shifted =[ data_shifted ; Ymuscles(:,:,i)'+c+r];  % Concatenating the data
    
    % Checking for colinearity and correlation between the regresssions
%     temp5=corrcoef(model{i}.C); % correlation coeficient
%     correlation(i,1)=temp5(2);
%     temp10(i,:)=vif([model{i}.C nanmean(Ymuscles(481:491,:,i))']); %Variance inflation
%     temp5=vif([model{i}.C]);
%     impact(i,:)=temp5;
    
  %fit linear model 
    mdl{i}= fitlm(unit(:,:,i),nanmean(Ymuscles(epochIdxToPlot,:,i))','VarNames',{'Reactive','Contextual', labels(i).Data(1:end-1)},'Intercept',false);
    mdl2{i}= fitlm(unit(:,:,i),nanmean(Ymuscles(epochIdxToPlot,:,i))','VarNames',{'Reactive','Contextual', labels(i).Data(1:end-1)},'Intercept',true);
    if showPlot==1
        subplot 211
        hold on
        li{subj}=scatter(i,squeeze(nanmean(temp2(epochIdxToPlot,1,i))),100,"filled",'MarkerFaceColor', colorOrder(subj,:));
        errorbar(i,mdl{i}.Coefficients.Estimate(1),mdl{i}.Coefficients.SE(1),'k')
        
        subplot 212
        hold on
        li2{subj}=scatter(i,squeeze(nanmean(temp2(epochIdxToPlot,2,i))),100,"filled",'MarkerFaceColor', colorOrder(subj,:));
        errorbar(i,mdl{i}.Coefficients.Estimate(2),mdl{i}.Coefficients.SE(2),'k')
    end
    reactive_trace(:,i,subj)= temp2(:,1,i);%strides x muscles x subjects, weights for reactive
    contextual_trace(:,i,subj)= temp2(:,2,i); %strides x muscles x subject, weights for contextual
end
mdl_perSubject(subj,:) = mdl; %subject x 28 muscles
C_perSubject(:,:,subj) = C_indv; %336 x 2C x subjects
yhat_perSubject(:,:,subj) = reconstruction_indv; %336 x strides x subject.
yhatContext_perSubject(:,:,subj) = yhatContext; %336 x strides
yhatReactive_perSubject(:,:,subj) = yhatReactive; %336 x strides
data_perSubject(:,:,subj) = data; %336xstrides

%shifted
yhat_perSubject_shifted(:,:,subj) = reconstruction_indv_shifted; %336 x strides x subject.
yhatContext_perSubject_shifted(:,:,subj) = yhat_context_shifted; %336 x strides
yhatReactive_perSubject_shifted(:,:,subj) = yhat_reactive_shifted; %336 x strides
data_perSubject_shifted(:,:,subj) = data_shifted; %336xstrides

% index=find(impact>5); % fidn the regressor with high colinearity
% 
% Uf=Uf(:,1:size(temp2,1));
% set(gcf,'color','w')
Uf=Uf(:,1:size(temp2,1));

       
        if subj==n_sub && showPlot==1
            subplot 211
            ylabel({'Reactive';'AU'})
            yline(0)
            set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
            legend([li{:}],subID, 'location','Best')
            set(gcf,'color','w')
            
            subplot 212
            ylabel({'Contextual';'AU'})
            yline(0)
            set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
            legend([li2{:}],subID, 'location','Best')
            set(gcf,'color','w')
            if adapt
                sgtitle('Early Adapt')
            else
                sgtitle('Early OGPost')
            end
        end

end
set(findall(gcf,'-property','FontSize'),'FontSize',15)

%% get the VAF per subject
R2 = struct();

for shifted = 0:1 %no shfit first, then shift
    if shifted
        yhatContextCurr = yhatContext_perSubject_shifted;
        dataCurr = data_perSubject_shifted;
        yhatCurr = yhat_perSubject_shifted;
        yhatReactiveCurr = yhatReactive_perSubject_shifted;
    else
        yhatContextCurr = yhatContext_perSubject;
        dataCurr = data_perSubject;
        yhatCurr = yhat_perSubject;
        yhatReactiveCurr = yhatReactive_perSubject;
    end
    yhatContextPost = squeeze(nanmean(yhatContextCurr(:,epochIdxToPlot,:),2));
    dataPost =  squeeze(nanmean(dataCurr(:,epochIdxToPlot,:),2)); %336 x subjects
    yhatPost = squeeze(nanmean(yhatCurr(:,epochIdxToPlot,:),2)); %336 x subjects
    yhatReactivePost = squeeze(nanmean(yhatReactiveCurr(:,epochIdxToPlot,:),2)); %336 x subjects
    
    for relativeToMean = 0:1 %0 center first, then mean center
        R2Curr = [];
        for i = 1:size(dataPost,2)
            R2Curr(1,i) = my_Rsquared_coeff(dataPost(:,i), yhatPost(:,i), relativeToMean);
            R2Curr(2,i) = my_Rsquared_coeff(dataPost(:,i), yhatContextPost(:,i), relativeToMean);
            R2Curr(3,i) = my_Rsquared_coeff(dataPost(:,i), yhatReactivePost(:,i), relativeToMean);
        end
        eval(['R2.shifted' num2str(shifted) '.relativeToMean', num2str(relativeToMean) '=R2Curr']);
    end
end

%% compute the norm of yhatReactive, yhatContext, normalized by y actual norm, flag when ynorm is smaller than yhatnorm
%the length of each element in the regressor is computed as a way to
%quantify the reactive vs contextual component/contribution in
%reconstructing the data.
yhatContextPost = squeeze(nanmean(yhatContext_perSubject(:,epochIdxToPlot,:),2));
dataPost =  squeeze(nanmean(data_perSubject(:,epochIdxToPlot,:),2)); %336 x subjects
yhatPost = squeeze(nanmean(yhat_perSubject(:,epochIdxToPlot,:),2)); %336 x subjects
yhatReactivePost = squeeze(nanmean(yhatReactive_perSubject(:,epochIdxToPlot,:),2)); %336 x subjects

regNorm.yhatContextPostNorm = vecnorm(yhatContextPost);
regNorm.yhatReactivePostNorm = vecnorm(yhatReactivePost);
regNorm.dataPostNorm = vecnorm(dataPost);
regNorm.yhatPostNorm = vecnorm(yhatPost);

nonNanNorm.yhatContextPostNorm = nanvecnorm(yhatContextPost);
nonNanNorm.yhatReactivePostNorm = nanvecnorm(yhatReactivePost);
nonNanNorm.dataPostNorm = nanvecnorm(dataPost);
nonNanNorm.yhatPostNorm = nanvecnorm(yhatPost);

regNorm.flagDataLessThanHat = regNorm.dataPostNorm < regNorm.yhatPostNorm;%||y|| < ||yhat||
regNorm.flagDesp = 'norm of actual y < norm yhat, suggesting non-additive regressors to get the yhat';
regNorm.contextRatio = regNorm.yhatContextPostNorm./regNorm.dataPostNorm;
regNorm.reactiveRatio = regNorm.yhatReactivePostNorm./regNorm.dataPostNorm;

nonNanNorm.flagDataLessThanHat = nonNanNorm.dataPostNorm < nonNanNorm.yhatPostNorm;%||y|| < ||yhat||
nonNanNorm.contextRatio = nonNanNorm.yhatContextPostNorm./nonNanNorm.dataPostNorm;
nonNanNorm.reactiveRatio = nonNanNorm.yhatReactivePostNorm./nonNanNorm.dataPostNorm;

if groupAnalysis
    savePrefix = 'Group';
else
    savePrefix = 'IndivSub';
end
if adapt
    savePrefix = [savePrefix '_eaAdapt'];
else
    savePrefix = [savePrefix '_eaOGPost'];
end
% save([groupID '_' savePrefix '_indivMuslce_OGPost_NormRatio.mat'],'yhatContextPost','yhatReactivePost','dataPost',...
%     'yhatPost','flaggedNorm','contextRatio','reactiveRatio','flagDesp')

%% save data
mdlFitRes = struct();
mdlFitRes.R2 = R2;
mdlFitRes.R2.labels = {'Total','Context','Reactive'};
mdlFitRes.yhat = yhat_perSubject;
mdlFitRes.yhatContext = yhatContext_perSubject;
mdlFitRes.yhatReactive = yhatReactive_perSubject;
mdlFitRes.y = data_perSubject;
mdlFitRes.C = C_perSubject;
mdlFitRes.reactive = reactive_trace;
mdlFitRes.context = contextual_trace;
mdlFitRes.norm.regNorm = regNorm;
mdlFitRes.norm.nonNanNorm = nonNanNorm;
mdlFitRes.postIdx = epochIdxToPlot;
mdlFitRes.groupID = groupID;
mdlFitRes.muscleLabel = ytl;
mdlFitRes.COrder = {'Reactive','Cotextual'};
mdlFitRes.mdl = mdl_perSubject;
% if groupAnalysis && adapt
%     groupAda = mdlFitRes;
%     clear indivSubAda
% elseif groupAnalysis && ~adapt
%     groupOGPost = mdlFitRes;
%     clear indivSubAda
% elseif ~groupAnalysis && ~adapt
%     indivSubOGPost = mdlFitRes;
%     clear indivSubAda
% end

save([groupID '_' savePrefix '_modelFitRes.mat'],'mdlFitRes')

% save('Group_indivMuslce_modelFitOGPost.mat.mat','groupOGPost')


%% plot to visualize the fit
subIdx = 1; %For BATS: 7 (low), 1 (middle), 11 (high R2); 
%for auf: 11 (hight), 18(low), 14(middle), 3 (lowish)

clims = [-1 1]*1;
f = figure('units','normalized','outerposition',[0 0 0.75 0.75]);
subplot(1,4,1)
% imagesc((squeeze(unit(:,1,:))'))
imagesc(reshape(C_perSubject(:,1,subIdx),12,28)',clims) 
title('Reactive')
fs = 10;
set(gca,'XTick',[],'YTick',yt,'YTickLabel',ytl,'FontSize',fs)
set(gca,'CLim',[-1 1]*1,'FontSize',fs); %making sure that the first plot color scheme goes [-1 1] and making the name of the labels larger

subplot(1,4,2)
% imagesc((squeeze(unit(:,2,:))'))
imagesc(reshape(C_perSubject(:,2,subIdx),12,28)',clims) 
title('Contextual')

subplot(1,4,3)
imagesc(reshape(dataPost(:,subIdx),[12,28])',clims) 
title('Y')

subplot(1,4,4)
imagesc(reshape(yhatPost(:,subIdx),[12,28])',clims)
title('Yhat')
set(gca,'YTickLabels',{}); %gca applies to the last subplot only

% %C before bias removal
% Coriginal = reshape(Cmuscles,2,[],1) + C(:,tmbase)';
% Coriginal = reshape(Coriginal,2,12,28);
% subplot(1,4,3)
% imagesc(squeeze(Coriginal(1,:,:))')
% title('ReactiveBiased')
% subplot(1,4,4)
% imagesc(squeeze(Coriginal(2,:,:))')
% title('ContextualBiased')

% subplot(1,4,3)
% imagesc(squeeze(nanmean(Yhat_TR(:,strideToPlot,:,iteIdx),2))') %avg over strides
% title('Yhat')
% 
% subplot(1,4,4)
% imagesc(squeeze(nanmean(Ymuscles_TR(strideToPlot,:,:),1))') %avg overstrides
% title('Y')

color4checkerboards

% Making the figure a bit prettier
fs=14; %font size
colormap(flipud(map)) %changing the color map to the one tha defined about, we flipup the matrix bc the code does L-R and we want R-L
% c=flipud('gray');
% colormap(flipud(gray))
set(gcf,'color','w'); %setting the background white, gcf applies to the whole figure
% set(gca,'CLim',[-1 1]*1)
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
colorbar %Showing the colormap bar

% compute R2 of the fit
sgtitle([groupID savePrefix ' Sub' num2str(subIdx) ' Fit R^2 = ' num2str(nanmean(R2.shifted0.relativeToMean1(1,subIdx)))])
exportgraphics(f,['/Volumes/Research/Shuqi/Grant/2023Summer/Results/fitCheckerboard/' groupID savePrefix 'Sub' num2str(subIdx) '.png'],'Resolution',300)

return
%%  Plotting 
model2{1}.C=C_indv; %Save regressors for model
model2{1}.Out= reconstruction_indv; %Save reconstruction of the indiivudal muscles to the model format for plotting 
analysis=0; %Flag asking if you want to run the analysis
isF=0; %Flag for fast leg


legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg2(model2{1},data,Uf,analysis,[],isF)
% legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model2{1},data,Uf,analysis,[],isF)
% save('BATS_indv_muscles.mat','X2asym')
% end

%% Plotting time course for individual muscles
analysis=0

for i=2
%pick the muscle that you want
    
    % Pick the data that you want to plot 
    %    Xasym=[temp2(1:40,:,i);nan(1,size(temp2,2));temp2(41:480,:,i);nan(1,size(temp2,2));temp2(481:end,:,i)];
    Xasym=[temp2(481:end,:,i)];
    %      Xasym=[temp2(1:200,:,i)];
    % Xasym=[X2asym(1:40,:);nan(1,size(X2asym,2));X2asym(41:80,:);...
    %     nan(1,size(X2asym,2));X2asym(81:520,:);nan(1,size(X2asym,2));X2asym(521:end,:)];
    % Xasym=[X2asym(681:843,:);nan(1,size(X2asym,2))];
    
    figure
    subplot(2,1,1)
    hold on
    % scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),10,'k','filled')
    scatter(1:length(movmean(Xasym(:,1),binwith)), movmean(Xasym(:,1),binwith),'filled','MarkerFaceColor',"#EDB120") %"#77AC30" )%
        
    % legend('Baseline','AutoUpdate','off')
    legend('Negative')
    title(labels(i).Data)
    %     legend('Removal Perturbation','AutoUpdate','off')
    % uistack(pp,'bottom')
    yline(0)
    ylabel({'Reactive';'(A.U)'})
    %     ylabel({'Removal';'(A.U)'})
    xlabel('strides')
    
    if size(temp2(:,:,i),2)>=2
        % figure
        subplot(2,1,2)
        hold on
        scatter(1:length(movmean(Xasym(:,2),binwith)), movmean(Xasym(:,2),binwith),'filled','MarkerFaceColor'," #00008B")
        %
        legend('Contextual')
        % % legend('Switch','AutoUpdate','off')
        % % uistack(pp,'bottom')
        yline(0)
        ylabel({'Contextual';'(A.U)'})
        xlabel('strides')
    end
    if i<=14
        isF=0;
    else
        isF=1;
    end
    set(gcf,'color','w')
    % Plot the time course plus the regressors 
    legacy_vizSingleModel_FreeModel_ShortAdaptation_IndvLeg(model{i},Ymuscles(:,:,i)',Uf,analysis,{labels(i).Data(2:end-1)},isF)
    set(gcf,'color','w')
    
    
end

%% Plotting all the muscles for individual participants 

muscleOrder={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT',...
    'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
ytl=([strcat('f',muscleOrder) strcat('s',muscleOrder)]);  %List of muscle
ytl(end:-1:1) = ytl(:);
yt=1:length(ytl);
fs=14;
c=colormap(jet(12));

%EMG Data 

if type==1 %Testing
    groupID='BATS'
elseif type==2 %Training
    groupID='BATR'
elseif type==3 %Both groups
     groupID='BAT'
end

files = dir ([groupID '*params.mat']);


n_subjects = size(files,1);

ii=0;
for i =1:n_subjects
    ii=1+ii;
    sub{ii} = files(i).name;
    subID{ii} = sub{ii}(1:end-10);
end

figure()

for s=1:12
    
    subplot 211
    hold on
    
    li{s}=scatter(1:28,reactive_trace(s,:),'filled','MarkerFaceColor', c(s,:));
    
    subplot 212
    hold on
    scatter(1:28,contextual_trace(s,:),'filled','MarkerFaceColor', c(s,:))

end

legend([li{:}],subID{:},'AutoUpdate','off')
subplot 211
hold on 
ylabel({'Reactive';'AU'})
yline(0)
set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)

subplot 212
hold on 
set(gca,'XTick',yt,'XTickLabel',ytl,'FontSize',10)
yline(0)
ylabel({'Contextual';'AU'})

set(gcf,'color','w')


%% Geting the average of the first 5 strides post-adaptation per muscle
figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(481:485,1,i)),nanmean(temp2(481:485,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})

xlim([-0.3 1.4])
ylim([-0.3 1.4])
set(gcf,'color','w')

%% Geeting the average of Early Adapt and lateAdapt
figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(41:45,1,i)),nanmean(temp2(41:45,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})
% axis square
xlim([-1 3])
ylim([-1 1])

title('EarlyAdapt')
set(gcf,'color','w')


figure
hold on
for i=1:28
    if i<15
        Li{1}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'b');
        text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
    else
        Li{2}=scatter(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),100,"filled",'MarkerFaceColor', 'r')  ;
        text(nanmean(temp2(440:480,1,i)),nanmean(temp2(440:480,2,i)),{labels(i).Data(2:end-1)})
    end
    
end
legend([Li{:}],['Slow';'Fast'])
ylabel({'Contextual';'A.U'})
xlabel({'Reactive';'A.U'})
axis square
title('LateAdapt')
set(gcf,'color','w')
xlim([-.1 .85])
ylim([-.1 .85])

%% Comparing individual analysis with the leg specific analysis

Cslow=[Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]./vecnorm([Casym(1:size(C,1)/2,1)  Casym(1:size(C,1)/2,2)]); %SLOW
Cfast=[Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]./vecnorm([Casym(1+size(C,1)/2:end,1)  Casym(1+size(C,1)/2:end,2)]); %FAST

Ymodel=Ymodel';
Yslow=Ymodel(1:680,1:size(Ymodel,2)/2); %SLOW
Yfast=Ymodel(1:680,size(Ymodel,2)/2+1:end); %FAST



Cinv=pinv(Cslow);
Xslow = Cinv*Yslow'; %x= y/C
hatYslow=  Cslow *Xslow ; %yhat = C

% hatYslow= hatYslow';

Cinv=pinv(Cfast);
Xfast = Cinv*Yfast'; %x= y/C
hatYfast=  Cfast *Xfast ; %yhat = C

% hatYfast=hatYfast;
Yfast=Yfast';
Yslow=Yslow';

% Variance explained

ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
mid=ones(1,3);
N=100;
gamma=1.5; %gamma > 1 expands the white (mid) part of the map, 'hiding' low values. Gamma<1 does the opposite
gamma=1;
map=[flipud(mid+ (ex1-mid).*([1:N]'/N).^gamma); mid; (mid+ (ex2-mid).*([1:N]'/N).^gamma)];



%% Variance all muscles
Rsquared_slow2= my_Rsquared_coeff(Yslow,hatYslow);
Rsquared_fast2 = my_Rsquared_coeff(Yfast,hatYfast);

Rsquared_slow_unit2 = my_Rsquared_coeff(Yslow,reconstruction_indv(1:168,:));
Rsquared_fast_unit2 = my_Rsquared_coeff(Yfast,reconstruction_indv(169:end,:));

figure()
showPlot(Rsquared_slow2)
hold on
showPlot(Rsquared_slow_unit2)
legend('Per leg','Per muscle')

figure
showPlot(Rsquared_fast2)
hold on
showPlot(Rsquared_fast_unit2)
legend('Per leg','Per muscle')

% both legs
Rsquared_both= my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],[hatYslow(:,:);hatYfast(:,:)]);
Rsquared_both_unit = my_Rsquared_coeff([Yslow(:,:);Yfast(:,:)],reconstruction_indv(:,:));

figure
showPlot(Rsquared_both)
hold on
showPlot(Rsquared_both_unit)
legend('Per leg','Per muscle')
title('Both leg')
set(gcf,'color','w')
% figure(3)
% ytl={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
%%
%% Calculating the R squared per muscle
% mm= 0:12:168;
% mm2=1:12:168;
% for m=1:14
%     Rsquared_slow(m,:) = my_Rsquared_coeff(Yslow(mm2(m):mm(m+1),:),hatYslow(mm2(m):mm(m+1),:));
%     Rsquared_fast(m,:) = my_Rsquared_coeff(Yfast(mm2(m):mm(m+1),:),hatYfast(mm2(m):mm(m+1),:));
%     
%     Rsquared_slow_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m)',Y2asym(:,:,m));
%     Rsquared_fast_unit(m,:) = my_Rsquared_coeff(Ymuscles(:,:,m+14)',Y2asym(:,:,m+14));
%     
% end
% % figure
% muscleOrder={'GLU','TFL','HIP','RF','VL','VM','BF', 'SEMB','SEMT','MG','LG','SOL','PER','TA'};
% ytl= defineMuscleListV2(muscleOrder); %List of muscle
% % ytl={'TA', 'PER', 'SOL', 'LG', 'MG', 'BF', 'SEMB', 'SEMT', 'VM', 'VL', 'RF', 'HIP','TFL', 'GLU'};
% % ytl(end:-1:1) = ytl(:);
% binw=5;
% 
% mtp=10;
% for m=mtp
%     %     figure
%     %      imagesc(nanmean(Yslow(mm2(m):mm(m+1),481:485),2)',[-1 ,1])
%     %      colormap(flipud(map))
%     figure
%     %     subplot(14,1,m)
%     hold on
%     aux1=conv(Rsquared_slow(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
%     
%     aux1=conv(Rsquared_slow_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
%     
%     %     plot(movmean(muscles_slow(m,:),5))
%     %     plot(movmean(muscles_slow_unit(m,:),5))
%     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
%     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
%     legend('Location','NorthEastOutside','AutoUpdate','off')
%     ylabel(ytl{m+14})
%     %     yline(nanmean(muscles_slow(m,10:30)))
%     
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     uistack(pp,'bottom')
% end
% set(gcf,'color','w')
% 
% % figure
% for m=mtp
%     %     subplot(14,1,m)
%     figure
%     hold on
%     aux1=conv(Rsquared_fast(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per leg','Color',"#0072BD") ;
%     
%     aux1=conv(Rsquared_fast_unit(m,:),ones(1,binw)/binw,'valid'); %Smoothing
%     plot(aux1,'LineWidth',2,'DisplayName','Per muscle','Color',"#A2142F") ;
%     %         plot(movmean(muscles_fast(m,:),5))
%     %         plot(movmean(muscles_fast_unit(m,:),5))
%     %     scatter(1:length(muscles(m,:)),movmean(muscles(m,:),5),'filled')
%     %     legend({'Per leg';'Per muscle'},'AutoUpdate','off');
%     
%     legend('Location','NorthEastOutside','AutoUpdate','off')
%     ylabel(ytl{m})
%     %     yline(nanmean(muscles_fast(m,10:30)))
%     pp=patch([40 480 480 40],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     %     pp=patch([0 440  440 0],[-0.5 -0.5 1 1],.7*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
%     
%     uistack(pp,'bottom')
% end
% set(gcf,'color','w')
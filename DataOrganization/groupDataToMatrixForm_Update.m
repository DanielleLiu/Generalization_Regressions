function [Y,Yasym,U,Ubreaks,Ysum,Yindv,labels,C,Cinv]=groupDataToMatrixForm_Update(subjIdx,fName,sqrtFlag)
% Original function created by Pablo Iturralde
% Update to work for the EMG generalization study. 
%Main update is adding the file name
%modified by DMMO 07/07/2023

%% Load real data:
EMGdata=h5read(fName,'/EMGdata');
SLA=h5read(fName,'/SLA');
speedDiff=h5read(fName,'/speedDiff');
breaks=h5read(fName,'/breaks');
labels=hdf5read(fName,'/labels');
regressors=h5read(fName,'/Regressors');

%%
U=speedDiff;
Ubreaks=breaks;

%%  Getting data 
muscPhaseIdx=1:size(EMGdata,2); %gettign the lenght of the muscle vector 
Y=EMGdata(:,muscPhaseIdx,subjIdx); %choosing the muscles and the participants that we want
Yindv=Y; %all data (steps x muslces X number of subjs)
Cinv=regressors(:,muscPhaseIdx,subjIdx);
if nargin>2 && sqrtFlag
    Y=sqrt(Y);
end
%Computes asymmety measure 
Yasym=Y-fftshift(Y,2);
Yasym=Yasym(:,1:size(Yasym,2)/2,:);

% Computer the sum the muscle activity per muscle
Ysum=Y+fftshift(Y,2);
Ysum=Ysum(:,1:size(Ysum,2)/2,:);

%Median across subjs
Y=nanmedian(Y,3);
C=nanmedian(Cinv,3)';
end

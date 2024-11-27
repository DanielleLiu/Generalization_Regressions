%% Free model - Linear regression - Indv Legs

fname='dynamicsData_BATR_subj_12_RemoveBadMuscles0_splits_0.h5'; %this is created from preProcessingLinearModel.m
load BATR_12_IndvLegsC13_ShortPertubations_RemovedBadMuscle_0.mat; 


EMGdata=h5read(fname,'/EMGdata');
binwith=5;
[Y,Yasym,Ycom,U,Ubreaks,Ysum]=groupDataToMatrixForm(1:size(EMGdata,3),0,fname);

Uf=[U;ones(size(U))];

Yf2= Ysum./2  + Yasym./2 ;
Ys2= Ysum./2 -  Yasym./2 ;

C=[C(:,5) C(:,6)]; 


C=[C(1:size(C,1)/2,:) C(size(C,1)/2+1:end,:)]; % ['s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}']

% Slow leg regression 
Cs=[C(:,1:4)];
Cinv=pinv(Cs)';
Ys=Y(:,1:size(Y,2)/2);
Yf=Y(:,size(Y,2)/2+1:end);

%Slow side 
% Cinv=pinv(C)';
Xs = Ys*Cinv; %x= y/C
Ys_hat= Cs * Xs' ; %yhat = C   

figure 
subplot(3,1,1)
plot( movmean(Xs(:,1),binwith))
hold on
plot(movmean(Xs(:,2),binwith))
plot(movmean(Xs(:,3),binwith))
plot(movmean(Xs(:,4),binwith))
% % plot(movmean(Xs(:,5),binwith))
% plot(movmean(Xs(:,6),binwith))
hold on
% title('C_2 = Late Adaptation')
title('S -Individual leg analysis')
set(gcf,'color','w')
yline(0) 
% legend('s_{TMbase}','s_{reactive}','s_{adaptive}')
% legend('s_{TMbase}','s_{reactive}','s_{adaptive}','f_{reactive}')
%  legend('s_{TMbase}','s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}','f_{TMbase}')
% legend('s_{reactive}','s_{TMbase}','f_{reactive}','f_{TMbase}')
legend('s_{reactive}','s_{adaptive}','f_{reactive}','f_{adaptive}')
% legend('s_{reactive}','s_{adaptive}','f_{reactive}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
%  legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')


 
subplot(3,1,2)
R2=1 - sum((Ys'-Ys_hat).^2)./sum((Ys'- mean(Ys2')).^2);
plot(movmean(R2,binwith))
ylabel('R^{2}')

subplot(3,1,3)
RMSE= sqrt(mean((Ys-Ys_hat').^2,2));
plot(movmean(RMSE,binwith))
ylabel('RMES')

model.C=Cs;
% legacy_vizSingleModelMLMC_FreeModel(model,Ys',Uf)
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Ys',Uf,0)

%% % Fast side regression
Cf=[ C(:,3:4) C(:,1:2)];
% Cf=[C(:,1:4)];
model.C=Cf;
Cinv=pinv(Cf)';
Xf = Yf*Cinv; %x= y/C
Yf_hat= Cf * Xf' ; %yhat = C   
figure 
subplot(3,1,1)
plot( movmean(Xf(:,1),binwith))
hold on
% title('C_1 = Early Adaptation')

% subplot(2,1,2)
plot(movmean(Xf(:,2),binwith))
plot(movmean(Xf(:,3),binwith))
plot(movmean(Xf(:,4),binwith))
% plot(movmean(Xf(:,5),binwith))
% plot(movmean(Xf(:,6),binwith))
hold on
% title('C_2 = Late Adaptation')
title('F -Individual leg analysis')
set(gcf,'color','w')
yline(0)
%  legend('s_{reactive}','s_{adaptive}','s_{TMbase}','f_{reactive}','f_{adaptive}','f_{TMbase}')
% legend('s_{reactive}','s_{TMbase}','f_{reactive}','f_{TMbase}')
legend('f_{reactive}','f_{adaptive}','s_{reactive}','s_{adaptive}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
% legend('s_{reactive}','s_{adaptive}','s_{earlypost}','f_{reactive}','f_{adaptive}','f_{earlypost}')
%  legend('s_{reactive}','s_{earlypost}','f_{reactive}','f_{earlypost}')
% legend('f_{TMbase}','f_{reactive}','f_{adaptive}')
% legend('s_{reactive}','f_{reactive}','f_{adaptive}')
subplot(3,1,2)
R2=1 - sum((Yf'-Yf_hat).^2)./sum((Yf'- mean(Yf2')).^2);
plot(movmean(R2,binwith))
ylabel('R^{2}')

subplot(3,1,3)
RMSE=[];
hold on 
RMSE= sqrt(mean((Yf-Yf_hat').^2,2));
plot(movmean(RMSE,binwith))
ylabel('RMES')

% model.C=C;
% legacy_vizSingleModelMLMC_FreeModel(model,Yf',Uf)
legacy_vizSingleModel_FreeModel_ShortAdaptation(model,Yf',Uf,0)

%% Combining left and right regressors

% C=[C1 C2 C4]; 
% Cs=[C(1:size(C,1)/2,:) C(size(C,1)/2+1:end,1)]; % ['s_{reactive}','s_{adaptive}','s_{earlypost}']
% Cinv=pinv(Cs)';
% Ys=Y(:,1:size(Y,2)/2);
% 
% 
% %Slow side 
% % Cinv=pinv(C)';
% Xs = Ys*Cinv; %x= y/C
% Ys_hat= Cs * Xs' ; %yhat = C   
% figure 
% subplot(2,1,1)
% plot( movmean(Xs(:,1),binwith))
% hold on
% % title('C_1 = Early Adaptation')
% 
% 
% 
% % subplot(2,1,2)
% plot(movmean(Xs(:,2),binwith))
% plot(movmean(Xs(:,3),binwith))
% plot(movmean(Xs(:,4),binwith))
% hold on
% % title('C_2 = Late Adaptation')
% title('S -Individual leg analysis')
% set(gcf,'color','w')
% yline(0)
% legend('s_{reactive}','s_{adaptive}','s_{EarlyPost}')
% legend('s_{reactive}','s_{adaptive}','s_{TMbase}','f_{reactive}')
% 
% 
% subplot(2,1,2)
% RMSE= sqrt(mean((Ys-Ys_hat').^2,2));
% 
% plot(movmean(RMSE,binwith))
% ylabel('RMES')
% 
% model.C=Cs;
% legacy_vizSingleModelMLMC_FreeModel(model,Ys',Uf)
% 
% 
% Cf=[C(size(C,1)/2+1:end,:) C(1:size(C,1)/2,1)];% ['f_{reactive}','f_{adaptive}','f_{earlypost}']
% Yf=Y(:,size(Y,2)/2+1:end);
% Cinv=pinv(Cf)';
% Xf = Yf*Cinv; %x= y/C
% Yf_hat= Cf * Xf' ; %yhat = C   
% figure 
% subplot(2,1,1)
% plot( movmean(Xf(:,1),binwith))
% hold on
% % title('C_1 = Early Adaptation')
% 
% % subplot(2,1,2)
% plot(movmean(Xf(:,2),binwith))
% plot(movmean(Xf(:,3),binwith))
% plot(movmean(Xf(:,4),binwith))
% hold on
% % title('C_2 = Late Adaptation')
% title('F -Individual leg analysis')
% set(gcf,'color','w')
% yline(0)
% legend('f_{reactive}','f_{adaptive}','f_{EarlyPost}')
% legend('f_{reactive}','f_{adaptive}','f_{TMbase}','s_{reactive}')
% 
% subplot(2,1,2)
% RMSE=[];
% RMSE= sqrt(mean((Yf-Yf_hat').^2,2));
% plot(movmean(RMSE,binwith))
% 
% ylabel('RMES')
% model.C=Cf;
% legacy_vizSingleModelMLMC_FreeModel(model,Yf',Uf)
% 

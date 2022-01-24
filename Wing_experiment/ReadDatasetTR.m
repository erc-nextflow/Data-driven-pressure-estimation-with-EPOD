%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPOD ESTIMATION OF FIELDS FROM WING EXPERIMENT %
% V. 1.000 - Stefano Discetti - 07/08/2021       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear, close all

fprintf('This code estimates velocity fields from the dataset of the wing \n')
fprintf('experiment performed by J.Chen and M. Raiola\n')

%---------------------------------------%
% Input data
%---------------------------------------%
Root='.\OUT_TRPIV\Wing_';  % Root of the field files (in .mat format)
NImg=3:10300;                                % Number of snapshots (in format [1:X]

step=30;                % downsampling step to build non-TR snapshots
TR=1:step:10000;         % Training snapshots (last one should be smaller than NImg(end)-tau-1
TE=1:10000;%1051:1250;           % Testing snapshots (last one should be smaller than NImg(end)-tau-1
tau=100;                % Length of the probe sequence
ROI=[41 110 41 150];    % Region of interest in the fields (number of vectors we are considering, indicated as first-last row, first-last column)
pry=50:5:100;           % position of the probes in the vertical direction according to the original PIV input file
prx=150;                % position of the probes in the horizontal direction according to the original PIV input file
Flagmean=1;             % if 1 it removes the mean in the data, if 0 it will not remove the mean
rmax=min([tau*numel(pry) numel(TR)]);   % max rank for reconstruction. Default option is min([tau*numel(pry) numel(TR)]);
% rmax = 50;
FlagFilt='';            % flag to activate the SG filter (active if 'YES')

prx=prx-ROI(3);
pry=pry-ROI(1);

% Reading input and computing mean (statistics only on Training dataset)
[U,V,PIV,PR,Hv,Wv,Nsnap]=ReadInput(Root,NImg,ROI,FlagFilt,TR,prx,pry,tau,Flagmean);

% Computing POD from training dataset
[PIV,PR]=CalcPOD(PIV,PR);

%Correlation matrix
Xi=PR.psi'*PIV.psi;
figure
imagesc(abs(Xi))
% Xi(abs(Xi) < (3/sqrt(length(Xi)))) = 0;
% Xi(abs(Xi) < 0.05) = 0;

%--------------------------%
% Testing                  %
%--------------------------%

% Read testing dataset
[PIV,PR]=CreateTestingDataset(U,V,TE,PIV,PR,Hv,Wv,tau,prx,pry,Flagmean);

PR.sigma(end,end)=PR.sigma(end-1,end-1)/100;
PIV.sigma(end,end)=PIV.sigma(end-1,end-1)/100;

% Temporal modes of the PROBE testing dataset
PR.psi_test=[PR.ute PR.vte]*PR.phi*PR.sigma^-1;

% Estimated temporal modes through EPOD
PIV.psi_est=(Xi\PR.psi_test')';
% PIV.psi_est=(Xi(1:rmax,1:rmax)\PR.psi_test(:,1:rmax)')';
% PIV.psi_est=(Xi(1:rmax,:)\PR.psi_test(:,1:rmax)')';

% Building reference velocity field and corresponding temporal modes
Uref=[reshape(PIV.ute',[numel(TE) Hv*Wv]) reshape(PIV.vte',[numel(TE) Hv*Wv])];
PIV.psi_ref=Uref*PIV.phi*PIV.sigma^-1;

% Comparison between rebuilt temporal modes, reference, and testing ones
figure
for i=1:16
    subplot(4,4,i)
        plot(TE,PIV.psi_ref(:,i),'r','Linewidth',1.5)
    hold on
        plot(TE,PIV.psi_est(:,i),'.b')
    ylim([-3/sqrt(numel(TR)) 3/sqrt(numel(TR))])
    plot(TR,PIV.psi(:,i),'s-k','Markersize',3)
    xlim([TE(1) TE(end)])
end

% Computing correlation coefficient between estimated and real time
% coefficients on testing samples
for i=1:size(PIV.psi_est,2)
    dum=corrcoef(PIV.psi_est(:,i),PIV.psi_ref(:,i));
    rho(i)=dum(1,2);
end
figure
plot(rho)
ylabel('$\rho$','interpreter','latex')


% Comparison between reconstructed velocity fields
[b,a] = butter(6,0.1); % parameter not adapted
% [b,a] = butter(6,0.5);
PIV.psi_est = filtfilt(b, a, PIV.psi_est);
% Urec=PIV.psi_est(:,1:rmax)*PIV.sigma(1:rmax,1:rmax)*PIV.phi(:,1:rmax)';
Urec=PIV.psi_est*PIV.sigma*PIV.phi';
UrefLO=PIV.psi_ref(:,1:rmax)*PIV.sigma(1:rmax,1:rmax)*PIV.phi(:,1:rmax)';

iFrame = 1450;

figure
subplot(3,2,1)
imagesc(reshape(Urec(iFrame,1:Hv*Wv),[Hv,Wv])+Flagmean*PIV.Um)
caxis([-6 6])
title('Reconstructed field - crosswise component')

subplot(3,2,2)
imagesc(reshape(Urec(iFrame,1+Hv*Wv:end),[Hv,Wv])+Flagmean*PIV.Vm)
caxis([0 18])

title('Reconstructed field - streamwise component')

subplot(3,2,3)
imagesc(reshape(UrefLO(iFrame,1:Hv*Wv),[Hv,Wv])+Flagmean*PIV.Um)
caxis([-6 6])
title('Low-order ground truth at rmax - crosswise component')

subplot(3,2,4)
imagesc(reshape(UrefLO(iFrame,1+Hv*Wv:end),[Hv,Wv])+Flagmean*PIV.Vm)
caxis([0 18])
title('Low-order ground truth at rmax - Streamwise component')

subplot(3,2,5)
imagesc(reshape(Uref(iFrame,1:Hv*Wv),[Hv,Wv])+Flagmean*PIV.Um)
caxis([-6 6])
title('Ground truth - crosswise component')

subplot(3,2,6)
imagesc(reshape(Uref(iFrame,1+Hv*Wv:end),[Hv,Wv])+Flagmean*PIV.Vm)
caxis([0 18])
title('Ground truth - streamwise component')

for i=1:6
    subplot(3,2,i)
    colormap jet(32)
    axis equal
end


% building error maps
err=sqrt(mean((Urec-Uref).^2,1));
errU=reshape(err(1:Hv*Wv),[Hv,Wv]);
errV=reshape(err(Hv*Wv+1:end),[Hv,Wv]);
stdref=max(std([PIV.utr' PIV.vtr'],1,1));


figure
subplot(1,2,1)
imagesc(errU./stdref)
colorbar
title('RMSE on U')
axis equal
colormap jet(32)
caxis([0 1])
subplot(1,2,2)
imagesc(errV./stdref)
colorbar
title('RMSE on V')
axis equal
colormap jet(32)
caxis([0 1])




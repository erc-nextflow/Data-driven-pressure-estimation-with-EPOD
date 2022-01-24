% edit on 20210513
% v11 trying more filters
% v12 adding gaussian noise to training data and probe data
% clc, clear, %close all
clear;
inputfolder='./DataInterp/';
outputfolder='./DataEstimation/';
if ispc
    inputfolder(inputfolder == '/') = '\';
    outputfolder(outputfolder == '/') = '\';
end
%%
Train=load([inputfolder,'TrainData']);
Probe=load([inputfolder,'ProbeData']);
X=Train.X;
Y=Train.Y;
% % adding noise
% Train.u = Train.u + 0.03*randn(size(Train.u));
% Train.v = Train.v + 0.03*randn(size(Train.v));
% Probe.u = Probe.u + 0.03*randn(size(Probe.u));
% Probe.v = Probe.v + 0.03*randn(size(Probe.v));

%%
Tau = 60;
if Train.Nimg(end)+Tau-1 > Probe.Nimg(end)
    disp('Train Set Exceed Probe Set Range from the End.');
    flag = find(Train.Nimg+Tau-1 < Probe.Nimg(end), 1, 'last');
    Train.Nimg = Train.Nimg(1:flag);
    Train.u = Train.u(:,1:flag);
    Train.v = Train.v(:,1:flag);
    Train.p = Train.p(:,1:flag);
end
clear Tau flag

%% Probe train snapshot matrix
Tau=60; %It's up to within how many frames will the flow pass by the field.
count=0; % for now count is i
clear ProbeTrain
for i=1:numel(Train.Nimg)
    count=count+1;
    for j=1:numel(Probe.Y)
        ind=find(Probe.Nimg==Train.Nimg(i));
        vec=(1+(j-1)*Tau):(j*Tau);
        ProbeTrain.u(vec,count)=Probe.u(j,ind+(0:Tau-1));
        ProbeTrain.v(vec,count)=Probe.v(j,ind+(0:Tau-1));
    end
    ProbeTrain.Nimg(count)=Train.Nimg(i);
end

%% Computation EPOD

UT=Train.u; VT=Train.v;
Um=mean(UT,2); UT=UT-Um;
Vm=mean(VT,2); VT=VT-Vm;
[PsiU,SigmaU,PhiU]=svd([UT; VT]','econ');

up=[ProbeTrain.u]; 
vp=[ProbeTrain.v];
um=mean(up,2); up=up-um;
vm=mean(vp,2); vp=vp-vm;
[PsiP,SigmaP,PhiP]=svd([up; vp]','econ');

Xi=PsiP'*PsiU;

figure; tmp = diag(SigmaU); plot(tmp(1:12)/sum(tmp), 'o-');
title('Singular Value in Flow Field Data');

figure; tmp = diag(SigmaP); plot(tmp(1:12)/sum(tmp), 'o-');
title('Singular Value in Probe Data');

% figure;
% subplot(2,3,1); pcolor(X, Y, reshape(PhiU(1:end/2,1),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; title('Mode 1');
% subplot(2,3,2); pcolor(X, Y, reshape(PhiU(1:end/2,2),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; title('Mode 2');
% subplot(2,3,3); pcolor(X, Y, reshape(PhiU(1:end/2,3),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; title('Mode 3');
% subplot(2,3,4); pcolor(X, Y, reshape(PhiU(1:end/2,4),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; title('Mode 4');
% subplot(2,3,5); pcolor(X, Y, reshape(PhiU(1:end/2,5),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; title('Mode 5');
% subplot(2,3,6); pcolor(X, Y, reshape(PhiU(1:end/2,6),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; title('Mode 6');
% sgtitle('First 6 POD Spatial Modes (u)');

figure; imagesc(abs(Xi(1:500, 1:500)));
axis equal
% colormap jet(16)

% figure; tmp = diag(SigmaU).^2; tmp = tmp/sum(tmp);
% plot(tmp(1:24), 'o-');
% tmp2 = 0*tmp;
% for itmp = 1:length(tmp)
%     tmp2(itmp) = sum(tmp(1:itmp));
% end
% hold on; yyaxis right;
% plot(tmp2(1:24), 'o-r');
% ylim([0 1])
% title('Singular Value in Flow Field Data');

% Xi(abs(Xi)<(3/(sqrt(rank(Xi)))))=0;
% Xi(abs(Xi)<(3/(sqrt(numel(Train.Nimg)))))=0;

%% Probe test snapshot matrix
Test=load([inputfolder,'TestData']);
count=0;
clear ProbeTest
for i=1:numel(Test.Nimg)-Tau
    count=count+1;
    for j=1:numel(Probe.Y)
        ind=find(Probe.Nimg==Test.Nimg(i));
        vec=(1+(j-1)*Tau):(j*Tau);
        ProbeTest.u(vec,count)=Probe.u(j,ind+(0:Tau-1));
        ProbeTest.v(vec,count)=Probe.v(j,ind+(0:Tau-1));
    end
    ProbeTest.Nimg(count)=Test.Nimg(i);
end

%% Computation of extended dynamics coefficients
psipr=[ProbeTest.u-mean(ProbeTest.u,2);ProbeTest.v-mean(ProbeTest.v,2)]'...
    *PhiP*(SigmaP^-1);
% sigma=sqrt(sum(psipr.^2,1));
% psipr=psipr/diag(sigma);
% psi_ext=(Xi\psipr')';
psi_ext=psipr*Xi;
% Nimg=ProbeTest.Nimg;
% save([outputfolder,'ExtendedPOD'],'Nimg','Tau','psi_ext','SigmaU',...
%     'PhiU','X','Y','Xi','Um','Vm','psipr');
% 

%% Plot extended temporal coefficients
% snap=98008;
snap = 96999;
ind=find(Test.Nimg==snap);

dt=0.1;
Dt=19; % depends on the sampling rate in real experiments
aTest=[Test.u-mean(Train.u,2);Test.v-mean(Train.v,2)]'*PhiU;
aExt=psi_ext*SigmaU;

nm=[1:4];
nt=(-80:1:80);
figure;
clf
plot(nt*dt, aTest(ind+nt,nm),'.')
hold on
plot(nt*dt,aExt(ind+nt,nm),'-')

xlabel('t')
ylabel('a')
% title('EPOD')

count=0;
indtrain=NaN*ones(size(1:Dt:numel(nt)));ttrain=[];
for i=1:Dt:numel(nt)
    count=count+1;
    try indtrain(count)=find(Train.Nimg==snap+nt(i)); end
    ttrain(count)=nt(i)*dt;
end
ttrain(isnan(indtrain))=[];
indtrain(isnan(indtrain))=[];
aTrain=PsiU(indtrain,:)*SigmaU;
plot(ttrain,aTrain(:,nm),'x');

legend('mode 1', 'mode 2', 'mode 3', 'mode 4', 'Location', 'southeast');

%% Display unit
iFrame = 98000;
flag1 = find(Test.Nimg == iFrame);
flag2 = find(Probe.Nimg== iFrame);
if isempty(flag1)|isempty(flag2)
    disp('No Relative Frame in Test Set or Probe Set');
elseif Probe.Nimg(end) < flag2+Tau-1
    disp('Exceed Probe Set');
else
    % edit one amount you want to display
    % Vdisp = sqrt(Test.u(:,flag1).^2+Test.v(:,flag1).^2);
    % Vdisp = abs(Test.u(:,flag1));
    Vdisp = Test.u(:,flag1);
    % Vdisp = Vdisp + 0.03*randn(size(Vdisp));
    % Vdisp = Test.u(:,fFrame1);
    Vdisp = reshape(Vdisp, size(Test.X));
    figure; pcolor(X, Y, Vdisp); shading interp; colormap(jet);
    axis equal; colorbar;
    title('Test set u');
    
    DispU = []; DispV = [];
    for j=1:numel(Probe.Y)
        DispU = [DispU; Probe.u(j,flag2+(0:Tau-1))'];
        DispV = [DispV; Probe.v(j,flag2+(0:Tau-1))'];
    end
    psiprII = [DispU-um;DispV-vm]'*PhiP/SigmaP;
    tmp = psiprII*Xi*SigmaU*PhiU';
    u_recon = tmp(1:length(tmp)/2) + Um';
    v_recon = tmp(length(tmp)/2+1:end) + Vm';
    u_recon = reshape(u_recon, size(Test.X));
    v_recon = reshape(v_recon, size(Test.X));
    % VdispII = sqrt(u_recon.^2+v_recon.^2);
    % VdispII = abs(u_recon);
    VdispII = u_recon;
    figure; pcolor(X, Y, VdispII); shading interp; colormap(jet);
    axis equal; colorbar;
    title('EPOD Estimation u');
    
    figure; pcolor(X, Y, VdispII - Vdisp); shading interp; colormap(jet);
    axis equal; colorbar;
    title('EPOD Estimation Error u');
    
end
clear u_recon v_recon Vdisp VdispII flag1 flag2 DispU DispV

%% Generating Estimated Data
FrameBegin = 95000; FrameEnd = 99800;
if FrameBegin < min(Probe.Nimg)
    disp('FrameBegin too Small');
end
if FrameEnd > max(Probe.Nimg)+Tau-1
    disp('FrameEnd too Big');
end
flag1 = find(FrameBegin <= Probe.Nimg, 1);
flag2 = find(FrameEnd >= Probe.Nimg, 1, 'last') + Tau - 1;

Nimg = Probe.Nimg(flag1:flag2-(Tau-1));
Recon_uv = [Probe.u(:,flag1:flag2); Probe.v(:,flag1:flag2)];
MidMatrix = PhiP/SigmaP*Xi*SigmaU;
BigMatrix = MidMatrix*PhiU';
save([outputfolder,'EPOD_Est.mat'],'Nimg','BigMatrix','um','vm','Um','Vm');
save([outputfolder,'EPOD_Est.mat'],'X','Y','Tau','Recon_uv','-append');
save([outputfolder,'EPOD_Est.mat'],'MidMatrix','PhiU','-append');
clear flag1 flag2

% post session of EPOD_from_probes for velocity and pressure estimation
% Woii junwei.chen@uc3m.es 210521
% v4 210607

clear
p_Wait = waitbar(0, 'Preparing...');
Est_path = 'EPOD_save.mat';
Save_path = 'EPOD_P.mat';
Test_Probe_path = 'Probes_Testing/Probes_';
Test_Field_path = 'Fields_Testing/Field_';
if ispc % change to Windows styled path
    Test_Probe_path(Test_Probe_path == '/') = '\';
    Test_Field_path(Test_Field_path == '/') = '\';
end
AFrame = 1:240;
Nu = 5e-5;
Rho = 1;
NoiseRatio = 0.00; % 0.03

load(Est_path);
Increx = X(1,2) - X(1,1);
Increy = Y(2,1) - Y(1,1);
Incret = 0.0065;
% reading probe data and reconstructing field
waitbar(0.05, p_Wait, 'Reading Probe Data for Testing...');
Upr = zeros(size(MidMatrix, 1)/2, length(AFrame));
Vpr = Upr;
iCount = 1;
for iFrame = AFrame
    load(sprintf('%s%06d.mat',Test_Probe_path,iFrame),'UP','VP');
    Upr(:, iCount) = UP(1:end/2)';
    Vpr(:, iCount) = VP(1:end/2)';
    iCount = iCount + 1;
end

% adding noise
Upr = Upr + NoiseRatio*randn(size(Upr));
Vpr = Vpr + NoiseRatio*randn(size(Vpr));

waitbar(0.15, p_Wait, 'Processing Probe Data for Testing...');
Um_pr = mean(Upr,2);
Vm_pr = mean(Vpr,2);
Upr = Upr - Um_pr;
Vpr = Vpr - Vm_pr;
% [psiP,sigmaP,phiP] = svd([Upr' Vpr'],'econ');
aExt = [Upr; Vpr]'*MidMatrix;
[b,a] = butter(6,0.1); % Butterworth filter on aExt, check it later
aExt = filtfilt(b, a, aExt);
F_est = aExt*phiF';
U_est = F_est(:, 1:end/2)'   + Um;
V_est = F_est(:,end/2+1:end)'+ Vm;
clear a b F_est Upr Vpr iCount iFrame MidMatrix UP VP

% reading test field data
U_ref = zeros(length(Um), length(AFrame));
V_ref = U_ref;
iCount = 1;
for iFrame = AFrame
    load(sprintf('%s%06d.mat',Test_Field_path,iFrame),'u','v');
    U_ref(:, iCount) = u(:);
    V_ref(:, iCount) = v(:);
    iCount = iCount + 1;
end

% estimating pressure
P_est = 0*U_est; P_ref = P_est;
parfor iFrame = AFrame(2):AFrame(end-1)
    u = reshape(U_ref(:,iFrame), size(X));
    v = reshape(V_ref(:,iFrame), size(X));
    ut= reshape(U_ref(:,iFrame+1)-U_ref(:,iFrame-1), size(X));
    vt= reshape(V_ref(:,iFrame+1)-V_ref(:,iFrame-1), size(X));
    p_init = PIterSolver(u,v,ut,vt,Nu,Rho,Increx,Increy,0*X);
    p_ref = p_init;
    P_ref(:,iFrame) = p_ref(:);
    u = reshape(U_est(:,iFrame), size(X));
    v = reshape(V_est(:,iFrame), size(X));
    ut= reshape(U_est(:,iFrame+1)-U_est(:,iFrame-1), size(X));
    vt= reshape(V_est(:,iFrame+1)-V_est(:,iFrame-1), size(X));
    p_est = PIterSolver(u,v,ut,vt,Nu,Rho,Increx,Increy,p_init);
    tmp1 = round(size(X,1)/4); tmp1 = tmp1:3*tmp1;
    tmp2 = round(size(X,2)/4); tmp2 = tmp2:3*tmp2;
    p_est=p_est-mean(p_est(tmp1,tmp2),'all')+mean(p_ref(tmp1,tmp2),'all');
    P_est(:,iFrame) = p_est(:);
end

close(p_Wait);

%% check time coefficient
figure;
hold on;
aTest = [U_ref-Um;V_ref-Vm]'*phiF;
plot(aTest(:,1),'-r')
plot(aExt(:,1),'-ob')

%% display single frame
iFrame = 12;
figure; pcolor(X, Y, reshape(U_est(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
title('U EPOD');
figure; pcolor(X, Y, reshape(U_ref(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
title('U REF');
figure; pcolor(X, Y, reshape(V_est(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
title('V EPOD');
figure; pcolor(X, Y, reshape(V_ref(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
title('V REF');

%% error map multi frames
U_err = reshape(sqrt(mean((U_est-U_ref).^2,2)),size(X));
V_err = reshape(sqrt(mean((V_est-V_ref).^2,2)),size(X));
figure; pcolor(X, Y, U_err);
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.3]);
title('U RMS ERROR');
figure; pcolor(X, Y, V_err);
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.15]);
title('V RMS ERROR');

%% display pressure in single frame
iFrame = 12;
figure; pcolor(X, Y, reshape(P_est(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar;
title('P EPOD');
figure; pcolor(X, Y, reshape(P_ref(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar;
title('P REF');
figure; pcolor(X, Y, reshape(P_est(:,iFrame)-P_ref(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar;
title('P ERROR');

%% display pressure in multi frame
P_err = P_est(:,2:end-1) - P_ref(:,2:end-1);
p_err = reshape(sqrt(mean(P_err.^2, 2)), size(X));
figure; pcolor(X, Y, p_err);
shading interp; colormap(jet); axis equal; colorbar; caxis([0, 0.02]);
title('P RMS ERROR');

%% functions
function D = LaplaceOperation2(f, dx, dy)
% 2D Laplacian operator on f
% calculates the values on the edges by linearly extrapolating the second
% differences from the interior
fx = 0.*f;
fx(:,2:end-1) = (f(:,1:end-2)+f(:,3:end)-2*f(:,2:end-1))./dx^2;
fx(:,1) = 2*fx(:,2) - fx(:,3);
fx(:,end) = 2*fx(:,end-1) - fx(:,end-2);
fy = 0.*f;
fy(2:end-1,:) = (f(1:end-2,:)+f(3:end,:)-2*f(2:end-1,:))./dy^2;
fy(1,:) = 2*fy(2,:) - fy(3,:);
fy(end,:) = 2*fy(end-1,:) - fy(end-2,:);
D = fx + fy;
end
function P = PIterSolver2(u,v,ut,vt,Nu,Rho,Increx,Increy,Pinit)
% Iterative Method for Pressure Estimation
% Pressure solver for time-resolved 2D2C velocity field
% need function LaplaceOperation2
% input
% u, v - velocity in current frame, arranged in meshgrid(X,Y)
% ut, vt - time derivative of u, v
% Nu, Rho - properties of fluid
% Increx, Increy - increment of meshgrid, uniform in the field
% Pinit - initial value of pressure, can be set to 0 in the first frame
[ux, uy] = gradient(u); ux = ux./Increx; uy = uy./Increy;
[vx, vy] = gradient(v); vx = vx./Increx; vy = vy./Increy;
delta_u = LaplaceOperation2(u, Increx, Increy);
delta_v = LaplaceOperation2(v, Increx, Increy);
px = -Rho.*(ut + u.*ux + v.*uy - Nu.*delta_u);
py = -Rho.*(vt + u.*vx + v.*vy - Nu.*delta_v);

P = zeros(size(u) + 2);
P(2:end-1, 2:end-1) = Pinit;
P(1,:) = P(2,:);   P(end,:) = P(end-1,:);
P(:,1) = P(:,2);   P(:,end) = P(:,end-1);
[xq, yq] = meshgrid([1.5:size(px,2)-0.5], 1:size(px,1));
px = interp2(px, xq, yq, 'spline');
px = [px(1,:);px;px(end,:)];
px1 = [0*px(:,1:2),px,0*px(:,1)];
px2 = [0*px(:,1),px,0*px(:,1:2)];
[xq, yq] = meshgrid(1:size(py,2), [1.5:size(py,1)-0.5]);
py = interp2(py, xq, yq, 'spline');
py = [py(:,1),py,py(:,end)];
py1 = [0*py(1:2,:);py;0*py(1,:)];
py2 = [0*py(1,:);py;0*py(1:2,:)];
Lambda = 0.1;
for iCount = 1:10000
    % Lambda = 0.2/sqrt(iCount) + 0.05; % Adaptive relaxation coefficient
    PD = Increx*px1 - (P - circshift(P, 1, 2))...
        -Increx*px2 - (P - circshift(P,-1, 2))...
        +Increy*py1 - (P - circshift(P, 1, 1))...
        -Increy*py2 - (P - circshift(P,-1, 1));
    P = P + Lambda.*PD;
    P(1,:) = P(2,:); P(end,:) = P(end-1,:);
    P(:,1) = P(:,2); P(:,end) = P(:,end-1);
    if mean(abs(PD),'all') < 0.0001
        break;
    end
end
P = P(2:end-1,2:end-1);
end

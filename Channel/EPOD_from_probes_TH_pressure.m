% post session of EPOD_from_probes_6 using Taylor's hypothesis to estimate
% pressure
% Woii junwei.chen@uc3m.es  210521
% v2 20210607 parellel computing pressure fix Incret

clear
Est_path = 'EPOD_save.mat';
Save_path = 'EPOD_P.mat';
Test_Field_path = 'Fields_Testing/Field_';
AFrame = 1:154;
Nu = 5e-5;
Rho = 1;

load(Est_path);
Increx = X(1,2) - X(1,1);
Increy = Y(2,1) - Y(1,1);
Incret = 0.0065;
clear a b F_est aExt Upr Vpr iCount iFrame MidMatrix phiF UP VP

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

% Taylor's hypothesis
U_tay = zeros(length(Um), length(AFrame));
V_tay = U_tay;
U_tay(:, 1) = U_ref(:, 1); V_tay(:, 1) = V_ref(:, 1);
iCount = 2;
u = reshape(U_tay(:, 1), size(X)); v = reshape(V_tay(:, 1), size(X));
for iFrame = AFrame(2:end)
%     % eular method
%     [ut, vt] = Taylor_Dt(u, v, X, Y);
%     u = u + ut*Incret;
%     v = v + vt*Incret;
%     U_tay(:, iCount) = u(:);
%     V_tay(:, iCount) = v(:);
%     iCount = iCount + 1;
    % Runge-Kutta 4th order
    [ut1, vt1] = Taylor_Dt(u, v, X, Y);
    [ut2, vt2] = Taylor_Dt(u+0.5*ut1*Incret, v+0.5*vt1*Incret, X, Y);
    [ut3, vt3] = Taylor_Dt(u+0.5*ut2*Incret, v+0.5*vt2*Incret, X, Y);
    [ut4, vt4] = Taylor_Dt(u+ut3*Incret, v+vt3*Incret, X, Y);
    u = u + (ut1+2*ut2+2*ut3+ut4)*Incret/6;
    v = v + (vt1+2*vt2+2*vt3+vt4)*Incret/6;
    U_tay(:, iCount) = u(:);
    V_tay(:, iCount) = v(:);
    iCount = iCount + 1;
end

% estimating pressure
% P_est = 0*U_est; P_ref = P_est;
P_ref = zeros(size(U_ref));
tmp1 = round(size(X,1)/4); tmp1 = tmp1:3*tmp1;
tmp2 = round(size(X,2)/4); tmp2 = tmp2:3*tmp2;
parfor iFrame = (AFrame(2):AFrame(end-1))-AFrame(1)+1
    u = reshape(U_ref(:,iFrame), size(X));
    v = reshape(V_ref(:,iFrame), size(X));
    ut= reshape(U_ref(:,iFrame+1)-U_ref(:,iFrame-1), size(X));
    vt= reshape(V_ref(:,iFrame+1)-V_ref(:,iFrame-1), size(X));
    p_init = PIterSolver(u,v,ut,vt,Nu,Rho,Increx,Increy,0*X);
    p_ref = p_init;
    P_ref(:,iFrame) = p_ref(:);
%     u = reshape(U_est(:,iFrame), size(X));
%     v = reshape(V_est(:,iFrame), size(X));
%     ut= reshape(U_est(:,iFrame+1)-U_est(:,iFrame-1), size(X));
%     vt= reshape(V_est(:,iFrame+1)-V_est(:,iFrame-1), size(X));
%     p_est = PIterSolver(u,v,ut,vt,Nu,Rho,Increx,Increy,p_init);
%     tmp1 = round(size(X,1)/4); tmp1 = tmp1:3*tmp1;
%     tmp2 = round(size(X,2)/4); tmp2 = tmp2:3*tmp2;
%     p_est=p_est-mean(p_est(tmp1,tmp2),'all')+mean(p_ref(tmp1,tmp2),'all');
%     P_est(:,iFrame) = p_est(:);
    u = reshape(U_tay(:,iFrame), size(X));
    v = reshape(V_tay(:,iFrame), size(X));
    [ut, vt] = Taylor_Dt(u, v, X, Y);
    p_tay = PIterSolver(u,v,ut,vt,Nu,Rho,Increx,Increy,p_init);
    p_tay=p_tay-mean(p_tay(tmp1,tmp2),'all')+mean(p_ref(tmp1,tmp2),'all');
    P_tay(:,iFrame) = p_tay(:);
end

%% display single frame
iFrame = 101;
% figure; pcolor(X, Y, reshape(U_est(:,iFrame),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
% title('U EPOD');
% figure; pcolor(X, Y, reshape(V_est(:,iFrame),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
% title('V EPOD');
figure; pcolor(X, Y, reshape(U_ref(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
title('U REF');
figure; pcolor(X, Y, reshape(V_ref(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
title('V REF');
figure; pcolor(X, Y, reshape(U_tay(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
title('U TAY');
figure; pcolor(X, Y, reshape(V_tay(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
title('V TAY');

%% error map multi frames
% U_err = reshape(sqrt(mean((U_est-U_ref).^2,2)),size(X));
% V_err = reshape(sqrt(mean((V_est-V_ref).^2,2)),size(X));
% figure; pcolor(X, Y, U_err);
% shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.3]);
% title('U RMS ERROR');
% figure; pcolor(X, Y, V_err);
% shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.15]);
% title('V RMS ERROR');

%% error map multi frames Taylor's hypo
U_err = reshape(sqrt(mean((U_tay-U_ref).^2,2)),size(X));
V_err = reshape(sqrt(mean((V_tay-V_ref).^2,2)),size(X));
figure; pcolor(X, Y, U_err);
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.3]);
title('U RMS ERROR TH');
figure; pcolor(X, Y, V_err);
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.15]);
title('V RMS ERROR TH');


%% display pressure in single frame
% iFrame = 101;
% figure; pcolor(X, Y, reshape(P_est(:,iFrame),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; caxis([-0.05 0.05]);
% title('P EPOD');
% figure; pcolor(X, Y, reshape(P_ref(:,iFrame),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; caxis([-0.05 0.05]);
% title('P REF');
% figure; pcolor(X, Y, reshape(P_est(:,iFrame)-P_ref(:,iFrame),size(X)));
% shading interp; colormap(jet); axis equal; colorbar; caxis([-0.05 0.05]);
% title('P ERROR');

%%
figure; pcolor(X, Y, reshape(P_tay(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.05 0.05]);
title('P TH');
figure; pcolor(X, Y, reshape(P_tay(:,iFrame)-P_ref(:,iFrame),size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.05 0.05]);
title('P ERROR TH');

%% display pressure in single frame
P_err = P_est(:,2:end-1) - P_ref(:,2:end-1);
p_err = reshape(sqrt(mean(P_err.^2, 2)), size(X));
figure; pcolor(X, Y, p_err);
shading interp; colormap(jet); axis equal; colorbar; caxis([0, 0.02]);
title('P RMS ERROR');

%%
P_err = P_tay(:,2:end-1) - P_ref(:,2:end-1);
p_err = reshape(sqrt(mean(P_err.^2, 2)), size(X));
figure; pcolor(X, Y, p_err);
shading interp; colormap(jet); axis equal; colorbar; caxis([0, 0.02]);
title('P RMS ERROR TH');

%% make video
v_target = VideoWriter('u_compare.avi');
v_target.FrameRate = 5;
open(v_target);
hf1 = figure; hold on;
for iFrame = 1:101
    subplot(1,3,1);     
    pcolor(X, Y, reshape(U_ref(:,iFrame),size(X)));
    shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
    xlim([0 1]); ylim([0 1]);
    set(gca, 'FontSize', 12);
    title('U REF', 'FontSize', 16);
    subplot(1,3,2);
    pcolor(X, Y, reshape(U_est(:,iFrame),size(X)));
    shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
    xlim([0 1]); ylim([0 1]);
    set(gca, 'FontSize', 12);
    title('U EPOD', 'FontSize', 16);
    subplot(1,3,3);
    pcolor(X, Y, reshape(U_tay(:,iFrame),size(X)));
    shading interp; colormap(jet); axis equal; colorbar; caxis([0 1.2]);
    xlim([0 1]); ylim([0 1]);
    set(gca, 'FontSize', 12);
    title('U TH', 'FontSize', 16);
    set(hf1, 'position', [10, 10, 1500, 500]);
    writeVideo(v_target, getframe(gcf));
end
close(v_target);
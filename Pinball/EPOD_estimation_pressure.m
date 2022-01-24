% post session of EPOD_estimation_12
% v8 branch: evaluate the accuracy between EPOD and TH
% Woii 210519
clear
Est_path = 'DataEstimation/EPOD_Est.mat';
Test_path = 'DataInterp/TestData.mat';
Save_path = 'DataEstimation/EPOD_P.mat';
if ispc
    Est_path(Est_path   == '/') = '\';
    Test_path(Test_path == '/') = '\';
    Save_path(Save_path == '/') = '\';
end
fBegin = 28;                   fEnd = 4525;
fIncre = 3;                    tIncre = 1;
Re = 130;
Nu = 1/Re;
Rho = 1;
flag_Iter = 1;
flag_nonEPOD = 1;
flag_SGF = 0; % Savitzky-Golay filter for probe data.
flag_Tfilt = 1; % filter to time coefficient
flag_SaveU = 1;
flag_Waitbar = 1;

if flag_Waitbar == 1
    p_Wait = waitbar(0, 'preparing...');
end
Est = load(Est_path);
if flag_Tfilt == 1
    aExt = EPODVeloFilGen(Est);
end
Test = load(Test_path);
if flag_SGF == 1
    Est.Recon_uv = sgolayfilt(Est.Recon_uv',2,5)';
end
X = Est.X;                     Y = Est.Y;
Increx = X(1,2) - X(1,1); % Increx = 0.08;
Increy = Y(2,1) - Y(1,1); % Increy = 0.08;
Incret = 0.1;
flag_FullRef = 1;
for iFrame = fBegin:fIncre:fEnd
    flag = find(Est.Nimg(iFrame) == Test.Nimg,1);
    if isempty(flag)
        disp('No Relevant Frame in TestData for the Reference');
        flag_FullRef = 0;
    end
end
if flag_Iter == 1
    p_EPOD_Iter = [];
    p_i = zeros(size(X));
end
if flag_SaveU == 1
    U_est = [];
    U_ref = [];
    Ut_est = [];
    Ut_ref = [];
end
flag_firstFrame = 1;
for iFrame = fBegin:fIncre:fEnd
    if flag_Waitbar == 1
        tmp1 = (fEnd-fBegin)/fIncre + 1;
        tmp2 = (iFrame-fBegin)/fIncre + 1;
        if flag_nonEPOD == 1
            waitbar(tmp2/tmp1/2, p_Wait,...
                ['estimating EPOD pressure fields... ',...
                num2str(tmp2),'/',num2str(tmp1)]);
        else
            waitbar(tmp2/tmp1, p_Wait,...
                ['estimating EPOD pressure fields... ',...
                num2str(tmp2),'/',num2str(tmp1)]);
        end
    end
    if flag_Tfilt == 1
        [u1, v1] = aExtVeloGen(aExt(iFrame-tIncre,:),...
            Est.PhiU,Est.Um,Est.Vm,Est.X); % previous  frame
        [u2, v2] = aExtVeloGen(aExt(iFrame,       :),...
            Est.PhiU,Est.Um,Est.Vm,Est.X); % current   frame
        [u3, v3] = aExtVeloGen(aExt(iFrame+tIncre,:),...
            Est.PhiU,Est.Um,Est.Vm,Est.X); % following frame
    else
        if fIncre == 1 & tIncre == 1
            if flag_firstFrame == 1
                [u1,v1] = EPODVeloGen(iFrame-tIncre,Est);% previous  frame
                [u2,v2] = EPODVeloGen(iFrame,       Est);% current   frame
                [u3,v3] = EPODVeloGen(iFrame+tIncre,Est);% following frame
                flag_firstFrame = 0;
            else
                u1 = u2; v1 = v2; u2 = u3; v2 = v3;
                [u3, v3] = EPODVeloGen(iFrame+tIncre,Est);
            end
        else
            [u1, v1] = EPODVeloGen(iFrame-tIncre,Est); % previous  frame
            [u2, v2] = EPODVeloGen(iFrame,       Est); % current   frame
            [u3, v3] = EPODVeloGen(iFrame+tIncre,Est); % following frame
        end
    end
    ut = (u3 - u1)./(2*Incret*tIncre);
    vt = (v3 - v1)./(2*Incret*tIncre);
    if flag_Iter == 1
        p_i = PIterSolver(u2, v2, ut, vt, Nu, Rho, Increx, Increy, p_i);
        if flag_FullRef == 1
            p_i = PConstCrt(p_i, X, Test, Est.Nimg(iFrame));
        end
        p_EPOD_Iter = [p_EPOD_Iter, p_i(:)];
    end
    if flag_SaveU == 1
        U_est = [U_est, [u2(:);v2(:)]];
        Ut_est = [Ut_est, [ut(:);vt(:)]];
    end
end

Nimg = Est.Nimg(fBegin:fIncre:fEnd);
p_Ref = [];
for iFrame = fBegin:fIncre:fEnd
    flag = find(Est.Nimg(iFrame) == Test.Nimg,1);
    p_Ref = [p_Ref, Test.p(:,flag)];
end
save(Save_path, 'Nimg', 'X', 'Y', 'p_Ref');

p_TH_Iter = [];
p_i = zeros(size(X));
for iFrame = fBegin:fIncre:fEnd
    tmp1 = (fEnd-fBegin)/fIncre + 1;
    tmp2 = (iFrame-fBegin)/fIncre + 1;
    waitbar(tmp2/tmp1, p_Wait,['estimating TH pressure fields... ',...
        num2str(tmp2),'/',num2str(tmp1)]);
    u2 = reshape(Test.u(:,iFrame), size(X));
    v2 = reshape(Test.v(:,iFrame), size(X));
%     % use mean velocity over frames
%     u_mean = reshape(Est.Um, size(X));
%     v_mean = reshape(Est.Vm, size(X));
    % use filtered current frame
    u_mean = imfilter(u2, fspecial('gaussian', [5, 9], 3), 'replicate');
    v_mean = imfilter(v2, fspecial('gaussian', [5, 9], 3), 'replicate');
    u_fluc = u2 - u_mean;
    v_fluc = v2 - v_mean;
    [upx, upy] = gradient(u_fluc); upx = upx./Increx; upy = upy./Increy;
    [vpx, vpy] = gradient(v_fluc); vpx = vpx./Increx; vpy = vpy./Increy;
    ut = - u_mean.*upx - v_mean.*upy;
    vt = - u_mean.*vpx - v_mean.*vpy;
    p_i = PIterSolver(u2, v2, ut, vt, Nu, Rho, Increx, Increy, p_i);
    
    p_ref = reshape(p_Ref(:, size(p_TH_Iter,2)+1), size(X));
    tmp1 = round(size(X,1)/4); tmp1 = tmp1:3*tmp1;
    tmp2 = round(size(X,2)/4); tmp2 = tmp2:3*tmp2;
    p_i = p_i - mean(p_i(tmp1,tmp2),'all') + mean(p_ref(tmp1,tmp2),'all');
    
    p_TH_Iter = [p_TH_Iter, p_i(:)];
end

if flag_nonEPOD == 1
    p_Iter = []; p_i = zeros(size(X));
    for iFrame = fBegin:fIncre:fEnd
        if flag_Waitbar == 1
            tmp1 = (fEnd-fBegin)/fIncre + 1;
            tmp2 = (iFrame-fBegin)/fIncre + 1;
            waitbar(tmp2/tmp1/2 + 0.5, p_Wait,...
                ['estimating reference pressure fields... ',...
                num2str(tmp2),'/',num2str(tmp1)]);
        end
        flag = find(Est.Nimg(iFrame-tIncre) == Test.Nimg,1);
        u1 = reshape(Test.u(:,flag), size(X));
        v1 = reshape(Test.v(:,flag), size(X));
        flag = find(Est.Nimg(iFrame+tIncre) == Test.Nimg,1);
        u3 = reshape(Test.u(:,flag), size(X));
        v3 = reshape(Test.v(:,flag), size(X));
        flag = find(Est.Nimg(iFrame) == Test.Nimg,1);
        u2 = reshape(Test.u(:,flag), size(X));
        v2 = reshape(Test.v(:,flag), size(X));
        ut = (u3 - u1)./(2*Incret*tIncre);
        vt = (v3 - v1)./(2*Incret*tIncre);
        
        p_i = PIterSolver(u2, v2, ut, vt, Nu, Rho, Increx, Increy, p_i);
        p_i = PConstCrt(p_i, X, Test, Est.Nimg(iFrame));
        p_Iter = [p_Iter, p_i(:)];
        if flag_SaveU == 1
            U_ref = [U_ref, [u2(:);v2(:)]];
            Ut_ref = [Ut_ref, [ut(:);vt(:)]];
        end
    end
    save(Save_path, 'p_Iter', '-append');
end
if flag_Iter == 1
    save(Save_path, 'p_EPOD_Iter', '-append');
end
if flag_SaveU == 1
    save(Save_path, 'U_est', 'U_ref', '-append');
end
if flag_Waitbar == 1
    close(p_Wait);
end

%% DISPLAY PART ONE FRAME
iFrame = 992;

p_Ref_iFrame = reshape(p_Ref(:,iFrame), size(X));
p_EPOD_Iter_iFrame = reshape(p_EPOD_Iter(:,iFrame), size(X));
% p_Iter_iFrame = reshape(p_Iter(:,iFrame), size(X));
p_TH_Iter_iFrame = reshape(p_TH_Iter(:,iFrame), size(X));

figure; pcolor(X, Y, p_Ref_iFrame);
shading interp; colormap(jet); axis equal; colorbar;
title('PRESSURE SIMULATION'); caxis([-0.5 0.5]);
figure; pcolor(X, Y, p_EPOD_Iter_iFrame);
shading interp; colormap(jet); axis equal; colorbar;
title('PRESSURE EPOD ITER'); caxis([-0.5 0.5]);
% figure; pcolor(X, Y, p_Iter_iFrame);
% shading interp; colormap(jet); axis equal; colorbar;
% title('PRESSURE ITER'); caxis([-0.5 0.5]);
figure; pcolor(X, Y, p_TH_Iter_iFrame);
shading interp; colormap(jet); axis equal; colorbar;
title('PRESSURE TH ITER'); caxis([-0.5 0.5]);

figure; pcolor(X, Y, p_EPOD_Iter_iFrame - p_Ref_iFrame);
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
title('PRESSURE ERROR EPOD ITER - SIMULATION');
% figure; pcolor(X, Y, p_EPOD_Iter_iFrame - p_Iter_iFrame);
% shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
% title('PRESSURE ERROR EPOD ITER - ITER');
% figure; pcolor(X, Y, p_Iter_iFrame - p_Ref_iFrame);
% shading interp; colormap(jet); axis equal; colorbar; caxis([-0.04 0.04]);
% title('PRESSURE ERROR ITER - SIMULATION');
figure; pcolor(X, Y, p_TH_Iter_iFrame - p_Ref_iFrame);
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2]);
title('PRESSURE ERROR TH ITER - SIMULATION');

%% MEAN ERROR MAP
figure; pcolor(X, Y, reshape(mean(abs(p_Iter-p_Ref),2), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.04]);
title('MEAN PRESSURE ERROR ITERATIVE - SIMULATION');
figure; pcolor(X, Y, reshape(mean(abs(p_EPOD_Iter-p_Ref),2), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.1]);
title('MEAN PRESSURE ERROR EPOD ITERATIVE - SIMULATION');
figure; pcolor(X, Y, reshape(mean(abs(p_EPOD_Iter-p_Iter),2), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.1]);
title('MEAN PRESSURE ERROR EPOD ITERATIVE - ITERATIVE');

%% MEAN ERROR MAP (RMS)
% figure; pcolor(X, Y, reshape(sqrt(mean((p_Iter-p_Ref).^2,2)), size(X)));
% shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.04]);
% title('MEAN PRESSURE ERROR ITERATIVE - SIMULATION');
figure;
pcolor(X, Y, reshape(sqrt(mean((p_EPOD_Iter-p_Ref).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.1]);
title('MEAN PRESSURE ERROR EPOD ITERATIVE - SIMULATION');
% figure;
% pcolor(X, Y, reshape(sqrt(mean((p_EPOD_Iter-p_Iter).^2,2)), size(X)));
% shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.1]);
% title('MEAN PRESSURE ERROR EPOD ITERATIVE - ITERATIVE');
figure;
pcolor(X, Y, reshape(sqrt(mean((p_TH_Iter-p_Ref).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.1]);
title('MEAN PRESSURE ERROR TH ITERATIVE - SIMULATION');

%% MEAN EPOD ERROR
if flag_SaveU ~= 1
    disp('TURN ON flag_SaveU');
end
figure;
pcolor(X, Y, reshape(sqrt(mean((U_est(1:end/2,:)...
    -U_ref(1:end/2,:)).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.25]);
title('MEAN EPOD ERROR u');
figure;
pcolor(X, Y, reshape(sqrt(mean((U_est(end/2+1:end,:)...
    -U_ref(end/2+1:end,:)).^2,2)), size(X))); caxis([0 0.25]);
shading interp; colormap(jet); axis equal; colorbar;
title('MEAN EPOD ERROR v');

figure;
pcolor(X, Y, reshape(sqrt(mean((Ut_est(1:end/2,:)...
    -Ut_ref(1:end/2,:)).^2,2)), size(X))); caxis([0 0.25]);
shading interp; colormap(jet); axis equal; colorbar;
title('MEAN EPOD ERROR \partial u/\partial t');
figure;
pcolor(X, Y, reshape(sqrt(mean((Ut_est(end/2+1:end,:)...
    -Ut_ref(end/2+1:end,:)).^2,2)), size(X))); caxis([0 0.25]);
shading interp; colormap(jet); axis equal; colorbar;
title('MEAN EPOD ERROR \partial v/\partial t');

[cu,cv,lu,lv] = SeePErrorComp(U_est,U_ref,X,Increx,Increy);
figure;pcolor(X, Y, reshape(sqrt(mean((cu).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.25]);
title('MEAN EPOD ERROR u convection');
figure;pcolor(X, Y, reshape(sqrt(mean((cv).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.25]);
title('MEAN EPOD ERROR v convection');
figure;pcolor(X, Y, Nu.*reshape(sqrt(mean((lu).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar;
title('MEAN EPOD ERROR u friction');
figure;pcolor(X, Y, Nu.*reshape(sqrt(mean((lv).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar;
title('MEAN EPOD ERROR v vfriction');

%% PRESSURE ERROR FRAME BY FRAME
v_target = VideoWriter('p_error.mp4');
v_target.FrameRate = 5;
open(v_target);
figure; set(gcf,'position',[100,100,1800,480]);
for iFrame = 1:size(p_Ref, 2)
    p_Ref_iFrame = reshape(p_Ref(:,iFrame), size(X));
    p_EPOD_Iter_iFrame = reshape(p_EPOD_Iter(:,iFrame), size(X));
    subplot(1,3,1);pcolor(X, Y, p_Ref_iFrame);
    shading interp; colormap(jet); axis equal; colorbar;
    title('PRESSURE SIMULATION'); caxis([-0.5 0.5]);
    subplot(1,3,2);pcolor(X, Y, p_EPOD_Iter_iFrame);
    shading interp; colormap(jet); axis equal; colorbar;
    title('PRESSURE EPOD ITERATIVE'); caxis([-0.5 0.5]);
    subplot(1,3,3);pcolor(X, Y, p_EPOD_Iter_iFrame - p_Ref_iFrame);
    shading interp; colormap(jet); axis equal; colorbar;
    title('PRESSURE ERROR'); caxis([-0.1 0.1]);
    set(gcf,'position',[100,100,1800,480]);
    writeVideo(v_target, getframe(gcf));
end
close(v_target);

%%
function [cu,cv,lu,lv] = SeePErrorComp(U_est,U_ref,X,Increx,Increy)
cu = []; cv = []; lu = []; lv = [];
for iFrame = 1:size(U_est, 2)
    u = reshape(U_est(1:end/2, iFrame), size(X));
    v = reshape(U_est(end/2+1:end, iFrame), size(X));
    [ux, uy] = gradient(u); ux = ux./Increx; uy = uy./Increy;
    [vx, vy] = gradient(v); vx = vx./Increx; vy = vy./Increy;
    delta_u = del2(u)*4./Increx; delta_v = del2(v)*4./Increx;
    cu_loop = u.*ux + v.*uy;
    cv_loop = u.*vx + v.*vy;
    
    u_ref = reshape(U_ref(1:end/2, iFrame), size(X));
    v_ref = reshape(U_ref(end/2+1:end, iFrame), size(X));
    [ux_ref, uy_ref] = gradient(u_ref);
    ux_ref = ux_ref./Increx; uy_ref = uy_ref./Increy;
    [vx_ref, vy_ref] = gradient(v_ref);
    vx_ref = vx_ref./Increx; vy_ref = vy_ref./Increy;
    delta_u_ref = del2(u_ref)*4./Increx;
    delta_v_ref = del2(v_ref)*4./Increx;
    cu_loop_ref = u_ref.*ux_ref + v_ref.*uy_ref;
    cv_loop_ref = u_ref.*vx_ref + v_ref.*vy_ref;
    
    cu = [cu, cu_loop_ref(:)-cu_loop(:)];
    cv = [cv, cv_loop_ref(:)-cv_loop(:)];
    lu = [lu, delta_u_ref(:)-delta_u(:)];
    lv = [lv, delta_v_ref(:)-delta_v(:)];
end
end
function p_i = PConstCrt(p_i, X, Test, FramePoint)
% substract the constant between p_i and refence
% nine points
tmp1 = floor(size(X,1)/4);
tmp2 = floor(size(X,2)/4);
tmp3 = mean(p_i(tmp1:tmp1:3*tmp1,tmp2:tmp2:3*tmp2), 'all');
flag = find(FramePoint == Test.Nimg,1);
p_refI = reshape(Test.p(:, flag), size(X));
tmp4 = mean(p_refI(tmp1:tmp1:3*tmp1,tmp2:tmp2:3*tmp2), 'all');
p_i = p_i - tmp3 + tmp4;
% five probes
% tmp1 = floor(size(X,1)/6);
% tmp3 = mean(p_i(tmp1:tmp1:5*tmp1,end), 'all');
% flag = find(FramePoint == Test.Nimg,1);
% p_refI = reshape(Test.p(:, flag), size(X));
% tmp4 = mean(p_refI(tmp1:tmp1:5*tmp1,end), 'all');
% p_i = p_i - tmp3 + tmp4;
% first point
% flag = find(FramePoint == Test.Nimg,1);
% p_refI = reshape(Test.p(:, flag), size(X));
% p_i = p_i - p_i(1) + p_refI(1);
end

function P = PIterSolver(u,v,ut,vt,Nu,Rho,Increx,Increy,Pinit)
% Iterative Method for Pressure Estimation
% Pressure solver for time-resolved 2D2C velocity field
% input
% u, v - velocity in current frame, arranged in meshgrid(X,Y)
% ut, vt - time derivative of u, v
% Nu, Rho - properties of fluid
% Increx, Increy - increment of meshgrid, uniform in the field
% Pinit - initial value of pressure, can be set to 0 in the first frame
[ux, uy] = gradient(u); ux = ux./Increx; uy = uy./Increy;
[vx, vy] = gradient(v); vx = vx./Increx; vy = vy./Increy;
delta_u = del2(u)*4./Increx^2; % to be edited when Increx ~= Increy
delta_v = del2(v)*4./Increx^2;
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
for iCount = 1:25000
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
% disp(iCount);
% disp(mean(abs(PD),'all'));
end

function [u, v] = EPODVeloGen(frameIndex, Est)
% Generating Velocity Field from EPOD
tmp = Est.Recon_uv(:, frameIndex:frameIndex+Est.Tau-1)';
Upr = (tmp(:)-[Est.um;Est.vm])';
Urecon = Upr*Est.BigMatrix + [Est.Um;Est.Vm]';
u = reshape(Urecon(1:end/2), size(Est.X));
v = reshape(Urecon(end/2+1:end), size(Est.X));
end

function aExt = EPODVeloFilGen(Est)
% Generating Velocity Field from EPOD, Filtering to aExt
[b,a] = butter(6,0.1);
nFrame = numel(Est.Nimg) - Est.Tau + 1;
% aExt = zeros(nFrame,size(Ext.MidMatrix,2));
aPr = zeros(size(Est.MidMatrix,1),nFrame);
for iFrame = 1:nFrame
    tmp = Est.Recon_uv(:, iFrame:iFrame+Est.Tau-1)';
    aPr(:,iFrame) = tmp(:);
    % aExt(iFrame,:) = (tmp(:)-[Est.um;Est.vm])'*Est.MidMatrix;
end
aExt = (aPr - [Est.um;Est.vm])'*Est.MidMatrix;
aExt = filtfilt(b, a, aExt);
end

function [u, v] = aExtVeloGen(aExt, PhiU, Um, Vm, X)
Urecon = aExt*PhiU' + [Um;Vm]';
u = reshape(Urecon(1:end/2), size(X));
v = reshape(Urecon(end/2+1:end), size(X));
end

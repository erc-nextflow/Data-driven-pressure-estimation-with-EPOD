clear
Test_path = 'DataInterp/TestData.mat';
if ispc
    Test_path(Test_path == '/') = '\';
end
fBegin = 28;                   fEnd = 4525;
fIncre = 3;
FArray = fBegin:fIncre:fEnd;
Re = 130;
Nu = 1/Re;
Rho = 1;    
Test = load(Test_path);
X = Test.X;                    Y = Test.Y;
Increx = X(1,2) - X(1,1); % Increx = 0.08;
Increy = Y(2,1) - Y(1,1); % Increy = 0.08;
Incret = 0.1;
NoiseLevel = 0.000;

% reference pressure (substract truncation error in grid convertion)
% can not cover full frames
p_ref = zeros(numel(X),numel(FArray));
% ut_Array = circshift(Test.u(:,3:end) - Test.u(:,1:end-2), 1, 2);
% vt_Array = circshift(Test.v(:,3:end) - Test.v(:,1:end-2), 1, 2);
tmp = zeros(numel(X),1);
ut_Array = [tmp, Test.u(:,3:end)-Test.u(:,1:end-2), tmp]./(2*Incret);
vt_Array = [tmp, Test.v(:,3:end)-Test.v(:,1:end-2), tmp]./(2*Incret);
parfor iFrame = 1:numel(FArray)
    u2 = reshape(Test.u(:,FArray(iFrame)), size(X));
    v2 = reshape(Test.v(:,FArray(iFrame)), size(X));
    ut = reshape(ut_Array(:,FArray(iFrame)), size(X));
    vt = reshape(vt_Array(:,FArray(iFrame)), size(X));
    p_i = PIterSolver(u2, v2, ut, vt, Nu, Rho, Increx, Increy, 0*u2);
    p_ref(:,iFrame) = p_i(:) - p_i(1) + Test.p(1,FArray(iFrame));
end

% add noise
Test.u = Test.u + NoiseLevel*randn(size(Test.u));
Test.v = Test.v + NoiseLevel*randn(size(Test.v));

p_tay = zeros(numel(X),numel(FArray));
parfor iFrame = 1:numel(FArray)
    u2 = reshape(Test.u(:,FArray(iFrame)), size(X));
    v2 = reshape(Test.v(:,FArray(iFrame)), size(X));
    [ut, vt] = Taylor_Dt(u2, v2, X, Y);
    p_i = PIterSolver(u2, v2, ut, vt, Nu, Rho, Increx, Increy, 0*u2);
    p_tay(:,iFrame) = p_i(:) - p_i(1) + Test.p(1,FArray(iFrame));
end

% Taylor's hypothesis propagation
THIncre = 60;
u_tay = Test.u; v_tay = Test.v;
for iFrame = fBegin:THIncre:fEnd
    uth = reshape(Test.u(:,iFrame), size(X));
    vth = reshape(Test.v(:,iFrame), size(X));
    for i_sub = 0:THIncre-1
        u_tay(:,iFrame+i_sub) = uth(:);
        v_tay(:,iFrame+i_sub) = vth(:);
        [ut1, vt1] = Taylor_Dt(uth, vth, X, Y);
        [ut2, vt2] = Taylor_Dt(uth+.5*ut1*Incret, vth+.5*vt1*Incret, X, Y);
        [ut3, vt3] = Taylor_Dt(uth+.5*ut2*Incret, vth+.5*vt2*Incret, X, Y);
        [ut4, vt4] = Taylor_Dt(uth+   ut3*Incret, vth+   vt3*Incret, X, Y);
        uth = uth + (ut1+2*ut2+2*ut3+ut4)*Incret/6;
        vth = vth + (vt1+2*vt2+2*vt3+vt4)*Incret/6;
    end
end

if NoiseLevel > 0
    u_tay(isnan(u_tay)) = 0; u_tay(isinf(u_tay)) = 0;
    v_tay(isnan(v_tay)) = 0; v_tay(isinf(v_tay)) = 0;
    % u_tay(isnan(u_tay)|isinf(u_tay)) = Test.u(isnan(u_tay)|isinf(u_tay));
    % v_tay(isnan(v_tay)|isinf(v_tay)) = Test.v(isnan(v_tay)|isinf(v_tay));
end
p_tay_adv = zeros(numel(X),numel(FArray));
parfor iFrame = 1:numel(FArray)
    u2 = reshape(u_tay(:,FArray(iFrame)), size(X));
    v2 = reshape(v_tay(:,FArray(iFrame)), size(X));
    [ut, vt] = Taylor_Dt(u2, v2, X, Y);
    % if NoiseLevel > 0
    %     ut(isnan(ut)) = 0; vt(isnan(vt)) = 0;
    % end
    p_i = PIterSolver(u2, v2, ut, vt, Nu, Rho, Increx, Increy, 0*u2);
    p_tay_adv(:,iFrame) = p_i(:) - p_i(1) + Test.p(1,FArray(iFrame));
end

%% DISPLAY PART ONE FRAME
iFrame = 992;
% p_ref_iFrame = reshape(Test.p(:,FArray(iFrame)), size(X));
p_ref_iFrame = reshape(p_ref(:,iFrame), size(X));
p_tay_iFrame = reshape(p_tay(:,iFrame), size(X));
p_tay_adv_iFrame = reshape(p_tay_adv(:,iFrame), size(X));

figure; pcolor(X, Y, p_ref_iFrame);
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.8 0.2])
AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('PRESSURE SIMULATION');

figure; pcolor(X, Y, p_tay_iFrame);
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.8 0.2])
AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('PRESSURE TH');

figure; pcolor(X, Y, p_tay_adv_iFrame);
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.8 0.2])
AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('PRESSURE TH ADVANCE');

figure; pcolor(X, Y, p_tay_iFrame - p_ref_iFrame);
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2])
AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('PRESSURE TH ERROR');

figure; pcolor(X, Y, p_tay_adv_iFrame - p_ref_iFrame);
shading interp; colormap(jet); axis equal; colorbar; caxis([-0.2 0.2])
AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('PRESSURE TH ADVANCE ERROR');

%% MEAN ERROR MAP
% p_ref = Test.p(:,FArray);

figure;
pcolor(X, Y, reshape(sqrt(mean((p_tay-p_ref).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.1]);
AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('MEAN PRESSURE ERROR TH');

figure;
pcolor(X, Y, reshape(sqrt(mean((p_tay_adv-p_ref).^2,2)), size(X)));
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.1]);
AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('MEAN PRESSURE ERROR TH ADVANCE');

%% ACCUMULATED ERROR ON TIME PROPAGATION
% ERR_TAY = sqrt(mean((p_tay - p_ref).^2, 1));
% ERR_TAY_ADV = sqrt(mean((p_tay_adv - p_ref).^2, 1));
% save('tmp.mat', 'ERR_TAY', 'ERR_TAY_ADV', '-append');

% ERR_EPOD = sqrt(mean((p_EPOD_Iter - p_ref).^2, 1));
% save('tmp.mat', 'ERR_EPOD', '-append');

t = 0:3:3*19;
figure; hold on;
semilogy(t, ERR_EPOD(1:20), 'o-');
semilogy(t, ERR_TAY(1:20), '^-');
semilogy(t, ERR_TAY_ADV(1:20), 'v-');
% semilogy(t, ERR_EPOD(1:20), 'o-', t, ERR_TAY(1:20), '^-',...
%     t, ERR_TAY_ADV(1:20), 'v-');
legend('Extended POD', 'Taylor''s hypothesis',...
    'Taylor''s hypothesis with propagation', 'Location', 'Northwest');
ylim([0 4])
xlabel('t - t_0'); ylabel('\epsilon');

%%
ERR_TAY = sqrt(mean((p_tay - p_ref).^2, 1));
ERR_TAY_ADV = sqrt(mean((p_tay_adv - p_ref).^2, 1));
t = 0:3:3*19;
ERR_KNN(1) = mean(ERR_KNN(2:end-1)); ERR_KNN(end) = mean(ERR_KNN(2:end-1));
ERR_KNN_TH(1) = mean(ERR_KNN_TH(2:end-1));
ERR_KNN_TH(end) = mean(ERR_KNN_TH(2:end-1));
figure; hold on;
plot(t, ERR_EPOD(1:20), 'o-');
plot(t, ERR_TAY(1:20), '^-');
plot(t, ERR_TAY_ADV(1:20), 'v-');
plot(t, ERR_KNN(1:20), '<-');
plot(t, ERR_KNN_TH(1:20), '>-');
% semilogy(t, ERR_EPOD(1:20), 'o-', t, ERR_TAY(1:20), '^-',...
%     t, ERR_TAY_ADV(1:20), 'v-');
legend('Extended POD', 'Taylor''s hypothesis',...
    'Taylor''s hypothesis with propagation', 'KNN',...
    'KNN with Taylor''s hypothesis correction', 'Location', 'Northwest');
ylim([1e-2 4])
xlabel('t - t_0'); ylabel(char(949));
set(gca,'Yscale','log')
set(gca,'YGrid','on')
set(gca, 'FontSize', 12);

%%
t = 0:3:3*19;
ERR_KNN(1) = mean(ERR_KNN(2:end-1)); ERR_KNN(end) = mean(ERR_KNN(2:end-1));
ERR_KNN_TH(1) = mean(ERR_KNN_TH(2:end-1));
ERR_KNN_TH(end) = mean(ERR_KNN_TH(2:end-1));
figure; hold on;
plot(t, ERR_EPOD(1:20), 'o-');
plot(t, ERR_KNN(1:20), '<-');
plot(t, ERR_KNN_TH(1:20), '>-');
% semilogy(t, ERR_EPOD(1:20), 'o-', t, ERR_TAY(1:20), '^-',...
%     t, ERR_TAY_ADV(1:20), 'v-');
legend('Extended POD', 'KNN',...
    'KNN with Taylor''s hypothesis correction', 'Location', 'Northwest');
ylim([1e-2 4])
xlabel('t - t_0'); ylabel(char(949));
set(gca,'Yscale','log')
set(gca,'YGrid','on')
set(gca, 'FontSize', 12);

%%
figure;
ERR_KNN(1) = mean(ERR_KNN(2:end-1)); ERR_KNN(end) = mean(ERR_KNN(2:end-1));
ERR_KNN_TH(1) = mean(ERR_KNN_TH(2:end-1));
ERR_KNN_TH(end) = mean(ERR_KNN_TH(2:end-1));
[mean(ERR_TAY_ADV) mean(ERR_TAY) mean(ERR_KNN) mean(ERR_KNN_TH) mean(ERR_EPOD)]
bar([mean(ERR_TAY_ADV) mean(ERR_TAY) mean(ERR_KNN) mean(ERR_KNN_TH) mean(ERR_EPOD)], 0.6, 'LineWidth', 0.1)
xlim([0.5 5.5]); ylim([0 0.16])
set(gca, 'FontSize', 12);
% set(gca, 'XTickLabel', {'Taylor''s hypothesis with propagation';'Taylor''s hypothesis';'';'';'Extended POD'})

%%
figure;
ERR_KNN(1) = mean(ERR_KNN(2:end-1)); ERR_KNN(end) = mean(ERR_KNN(2:end-1));
ERR_KNN_TH(1) = mean(ERR_KNN_TH(2:end-1));
ERR_KNN_TH(end) = mean(ERR_KNN_TH(2:end-1));
[mean(ERR_TAY_ADV) mean(ERR_TAY) mean(ERR_KNN) mean(ERR_KNN_TH) mean(ERR_EPOD)]
bar([mean(ERR_KNN) mean(ERR_KNN_TH) mean(ERR_EPOD)], 0.6, 'LineWidth', 0.1)
xlim([0.5 3.5]); ylim([0 0.08])
set(gca, 'FontSize', 12);
% set(gca, 'XTickLabel', {'Taylor''s hypothesis with propagation';'Taylor''s hypothesis';'';'';'Extended POD'})

%%
Test_path = 'DataInterp/TestData.mat';
if ispc
    Test_path(Test_path == '/') = '\';
end
fBegin = 28;                   fEnd = 4525;
fIncre = 3;
FArray = fBegin:fIncre:fEnd;
Test = load(Test_path);
X = Test.X;                    Y = Test.Y;
Increx = X(1,2) - X(1,1); % Increx = 0.08;
Increy = Y(2,1) - Y(1,1); % Increy = 0.08;
Incret = 0.08;
utR = 0*X;
for iFrame = 1:numel(FArray)
    u1 = reshape(Test.u(:,FArray(iFrame)-1), size(X));
    v1 = reshape(Test.v(:,FArray(iFrame)-1), size(X));
    u2 = reshape(Test.u(:,FArray(iFrame)), size(X));
    v2 = reshape(Test.v(:,FArray(iFrame)), size(X));
    u3 = reshape(Test.u(:,FArray(iFrame)+1), size(X));
    v3 = reshape(Test.v(:,FArray(iFrame)+1), size(X));
    [utT, vtT] = Taylor_Dt(u2, v2, X, Y);
    utS = (u3-u1)./(2*Incret);
    utR = (utT-utS).^2 + utR;
end
utR = sqrt(utR/numel(FArray));

figure;
pcolor(X, Y, utR);
shading interp; colormap(jet); axis equal; colorbar; caxis([0 0.25]);
% AddMark; xlabel('x'); ylabel('y'); set(gca, 'FontSize', 12);
title('\partial u/\partial t TH');

%% functions
function [ut, vt] = Taylor_Dt(u, v, x, y)
Increx = x(1,2) - x(1,1);
Increy = y(2,1) - y(1,1);
u_mean = imgaussfilt(u, 7);
v_mean = imgaussfilt(v, 7);
u_fluc = u - u_mean;
v_fluc = v - v_mean;
[ux, uy] = gradient(u_fluc); ux = ux./Increx; uy = uy./Increy;
[vx, vy] = gradient(v_fluc); vx = vx./Increx; vy = vy./Increy;
ut = -u_mean.*ux - v_mean.*uy;
vt = -u_mean.*vx - v_mean.*vy;
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
for iCount = 1:10000
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
function AddMark
hold on;
xlim([-2, 7]); ylim([-3 3]);
theta = 2*pi*(0:100)/100;
xCyl1 = -1.299 + 0.5*cos(theta);
yCyl1 = 0.5*sin(theta);
xCyl2 = 0.5*cos(theta); xCyl3 = xCyl2;
yCyl2 = 0.75 + 0.5*sin(theta);
yCyl3 =-0.75 + 0.5*sin(theta);
fill(xCyl1, yCyl1, 'w'); fill(xCyl2, yCyl2, 'w'); fill(xCyl3, yCyl3, 'w');
xProbe = [7 7 7 7 7]; yProbe = -2:1:2;
% plot(xProbe, yProbe, 'ok', 'MarkerEdgeColor', 'k',...
%     'MarkerFaceColor','w', 'MarkerSize', 6);
hold off;
end

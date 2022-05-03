% script to estimate pressure and propagate velocity field using TH
% Woii junwei.chen@alumnos.uc3m.es 20210922 v1.0
% Experimental Aerodynamics and Propulsion Lab
% of Univeridad Carlos III de Madrid
% post processing of ReadDatasetTR.m
% run ReadDatasetTR.m first then this, without clearing working space
% replacing for with parfor in line 78 to accelerate pressure calculation

% U[P]ref: velocity[pressure] from PIV
% U[P]rec: velocity[pressure] from EPOD
% U[P]pro: velocity[pressure] from Taylor's hypothesis propagation
% Ptay: pressure from PIV snapshot using Taylor's propogation

AFrame = 6000:6030;             % frame number in Uref and Urec interested
Nu = 1e-6;                      % kinetic viscosity coefficient of water
Rho = 1e3;                      % density of water
Mm_per_px_ratio = 0.1197;       % spatial resolution of photos
Sample_rate = 30;               % sample rate of PIV
Vector_Spacing = 10;            % vector per pixel in PIV
% increment in space and time
Increx = Vector_Spacing*Mm_per_px_ratio*1e-3;
Increy = Vector_Spacing*Mm_per_px_ratio*1e-3;
Incret = 1/Sample_rate;
% generating grid (unit: mm)
[X, Y] = meshgrid(0:ROI(4)-ROI(3),0:ROI(2)-ROI(1));
X = X*Vector_Spacing*Mm_per_px_ratio;
Y = Y*Vector_Spacing*Mm_per_px_ratio;
% velocity smoothing using 3D Gaussian smoothing
Uref1 = reshape(Uref(:,1:end/2), [size(Uref,1),size(X,1),size(X,2)]);
Uref1 = imgaussfilt3(Uref1, 0.75);
Uref2 = reshape(Uref(:,end/2+1:end), [size(Uref,1),size(X,1),size(X,2)]);
Uref2 = imgaussfilt3(Uref2, 0.75);
UrefN = [reshape(Uref1, [size(Uref,1),size(Uref,2)/2]),...
    reshape(Uref2, [size(Uref,1),size(Uref,2)/2])];
Urec1 = reshape(Urec(:,1:end/2), [size(Urec,1),size(X,1),size(X,2)]);
Urec1 = imgaussfilt3(Urec1, 1);
Urec2 = reshape(Urec(:,end/2+1:end), [size(Urec,1),size(X,1),size(X,2)]);
Urec2 = imgaussfilt3(Urec2, 1);
UrecN = [reshape(Urec1, [size(Urec,1),size(Urec,2)/2]),...
    reshape(Urec2, [size(Urec,1),size(Urec,2)/2])];
clear Uref1 Uref2 Urec1 Urec2

% velocity field propogation using Taylor's hypothesis
Upro = zeros(AFrame(end)-AFrame(1)+1,size(X,1),size(X,2));
Vpro = zeros(AFrame(end)-AFrame(1)+1,size(X,1),size(X,2));
u = reshape(UrefN(AFrame(1),end/2+1:end), size(X));
v = reshape(UrefN(AFrame(1),1:end/2), size(X));
u = (u + Flagmean*PIV.Vm).*Sample_rate.*Mm_per_px_ratio.*1e-3;
v = (v + Flagmean*PIV.Um).*Sample_rate.*Mm_per_px_ratio.*1e-3;
Upro(1,:,:) = u;
Vpro(1,:,:) = v;
iCount = 1;
for iFrame = 2:AFrame(end)-AFrame(1)+1
    iCount = iCount + 1;
    % unidirectional propagation using 4th-order Runge-Kutta method
    [ut1, vt1] = Taylor_Dt(u, v, X.*1e-3, Y.*1e-3);
    [ut2, vt2] = Taylor_Dt(u+0.5*ut1*Incret, v+0.5*vt1*Incret,...
        X.*1e-3, Y.*1e-3);
    [ut3, vt3] = Taylor_Dt(u+0.5*ut2*Incret, v+0.5*vt2*Incret,...
        X.*1e-3, Y.*1e-3);
    [ut4, vt4] = Taylor_Dt(u+ut3*Incret, v+vt3*Incret, X.*1e-3, Y.*1e-3);
    u = u + (ut1+2*ut2+2*ut3+ut4)*Incret/6;
    v = v + (vt1+2*vt2+2*vt3+vt4)*Incret/6;
    Upro(iCount,:,:) = u;
    Vpro(iCount,:,:) = v;
end
clear u v iCount iFrame ut1 ut2 ut3 ut4 vt1 vt2 vt3 vt4

if AFrame(1)-1 < 1
    fprintf('error:the frame previous to AFrame(1) should be in Uref\n');
end
if AFrame(end)+1 > size(Uref,1)
    fprintf('error:the frame next to AFrame(end) should be in Uref\n');
end
% Pressure estimation
Pref = zeros(numel(AFrame), size(Uref,2)/2);
Ptay = Pref; Prec = Pref; Ppro = Pref;
parfor iFrame = 1:length(AFrame)
    % u[v]1: previous frame, u[v]2: current frame, u[v]3: next frame;
    % u[v]t: time-derivative of velocity
    % Reference Pressure
    u1 = reshape(UrefN(AFrame(iFrame)-1,end/2+1:end), size(X));
    v1 = reshape(UrefN(AFrame(iFrame)-1,1:end/2), size(X));
    u2 = reshape(UrefN(AFrame(iFrame),end/2+1:end), size(X));
    v2 = reshape(UrefN(AFrame(iFrame),1:end/2), size(X));
    u3 = reshape(UrefN(AFrame(iFrame)+1,end/2+1:end), size(X));
    v3 = reshape(UrefN(AFrame(iFrame)+1,1:end/2), size(X));
    u1 = (u1 + Flagmean*PIV.Vm).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    u2 = (u2 + Flagmean*PIV.Vm).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    u3 = (u3 + Flagmean*PIV.Vm).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    v1 = (v1 + Flagmean*PIV.Um).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    v2 = (v2 + Flagmean*PIV.Um).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    v3 = (v3 + Flagmean*PIV.Um).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    ut = (u3 - u1)./(2*Incret);
    vt = (v3 - v1)./(2*Incret);
    p = PIterSolver2(u2, v2, ut, vt, Nu, Rho, Increx, Increy, 0*u2);
    Pref(iFrame, :) = p(:);
    % Taylor's hypothesis
    [ut, vt] = Taylor_Dt(u2, v2, X.*1e-3, Y.*1e-3);
    p = PIterSolver2(u2, v2, ut, vt, Nu, Rho, Increx, Increy, p);
    Ptay(iFrame, :) = p(:);
    % Pressure from EPOD
    u1 = reshape(UrecN(AFrame(iFrame)-1,end/2+1:end), size(X));
    v1 = reshape(UrecN(AFrame(iFrame)-1,1:end/2), size(X));
    u2 = reshape(UrecN(AFrame(iFrame),end/2+1:end), size(X));
    v2 = reshape(UrecN(AFrame(iFrame),1:end/2), size(X));
    u3 = reshape(UrecN(AFrame(iFrame)+1,end/2+1:end), size(X));
    v3 = reshape(UrecN(AFrame(iFrame)+1,1:end/2), size(X));
    u1 = (u1 + Flagmean*PIV.Vm).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    u2 = (u2 + Flagmean*PIV.Vm).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    u3 = (u3 + Flagmean*PIV.Vm).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    v1 = (v1 + Flagmean*PIV.Um).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    v2 = (v2 + Flagmean*PIV.Um).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    v3 = (v3 + Flagmean*PIV.Um).*Sample_rate.*Mm_per_px_ratio.*1e-3;
    ut = (u3 - u1)./(2*Incret);
    vt = (v3 - v1)./(2*Incret);
    p = PIterSolver2(u2, v2, ut, vt, Nu, Rho, Increx, Increy, p);
    Prec(iFrame, :) = p(:);
    % Pressure from TH propogation
    iFrame_tmp = (iFrame-1)*(AFrame(2)-AFrame(1))+1;
    u = squeeze(Upro(iFrame_tmp,:,:));
    v = squeeze(Vpro(iFrame_tmp,:,:));
    [ut, vt] = Taylor_Dt(u, v, X.*1e-3, Y.*1e-3);
    p = PIterSolver2(u, v, ut, vt, Nu, Rho, Increx, Increy, p);
    Ppro(iFrame, :) = p(:);
end
clear u1 u2 u3 v1 v2 v3 ut vt p iFrame iFrame_tmp

% interframe correction
[xcor0, ycor0] = meshgrid(1:size(X,2),1:size(X,1));
for iFrame = 2:length(AFrame)
    % fixed point
    % Pref(iFrame, :) = Pref(iFrame, :) - Pref(iFrame, 1) + Pref(1, 1);
    % minimum change of fluid parcels in adjacent frames
    u = Flagmean*PIV.Vm +...
        reshape(UrefN(AFrame(iFrame-1),end/2+1:end), size(X));
    v = Flagmean*PIV.Um +...
        reshape(UrefN(AFrame(iFrame-1),1:end/2), size(X));
    xcor1 = xcor0 - u./2.*(AFrame(2)-AFrame(1))./Vector_Spacing;
    ycor1 = ycor0 - v./2.*(AFrame(2)-AFrame(1))./Vector_Spacing;
    u = Flagmean*PIV.Vm +...
        reshape(UrefN(AFrame(iFrame),end/2+1:end), size(X));
    v = Flagmean*PIV.Um +...
        reshape(UrefN(AFrame(iFrame),1:end/2), size(X));
    xcor2 = xcor0 + u./2.*(AFrame(2)-AFrame(1))./Vector_Spacing;
    ycor2 = ycor0 + v./2.*(AFrame(2)-AFrame(1))./Vector_Spacing;
    p1 = interp2(reshape(Pref(iFrame-1,:),size(xcor1)), xcor1, ycor1);
    p2 = interp2(reshape(Pref(iFrame,:),size(xcor2)), xcor2, ycor2);
    Pref(iFrame,:) = Pref(iFrame,:) - mean(p2(:)-p1(:),'all','omitnan');
end
% The mean pressure in the domian should be equal to that from PIV velocity
% in all frames
for iFrame = 1:length(AFrame)
    Prec(iFrame, :) = Prec(iFrame, :)...
        - mean(Prec(iFrame, :)) + mean(Pref(iFrame,:));
    Ptay(iFrame, :) = Ptay(iFrame, :)...
        - mean(Ptay(iFrame, :)) + mean(Pref(iFrame,:));
    Ppro(iFrame, :) = Ppro(iFrame, :)...
        - mean(Ppro(iFrame, :)) + mean(Pref(iFrame,:));
end
clear u v xcor1 ycor1 xcor2 ycor2 p1 p2 minmin maxmax minmax
clear bregion bline1 bline2 sigma ipx bpeak1 bpeak2

%% show

% plot
% t = (AFrame - AFrame(1))./Sample_rate;
% figure; hold on;
% plot(t, sqrt(mean((Prec-Pref).^2,2)), '-o');
% plot(t, sqrt(mean((Ptay-Pref).^2,2)), '-^');
% plot(t, sqrt(mean((Ppro-Pref).^2,2)), '-v');
% legend('EPOD', 'Taylor', 'Taylor from t=0', 'Location', 'Northwest');
% 
% figure; pcolor(X, Y, reshape(Pref(1,:), size(X)))
% shading interp; colormap(jet); axis equal; colorbar; caxis([-1 1]);
% title('REF');
% figure; pcolor(X, Y, reshape(Prec(1,:), size(X)))
% shading interp; colormap(jet); axis equal; colorbar; caxis([-1 1]);
% title('EPOD');

grid = load('OUT_TRPIV/Grid_Wing.mat');

NACA0018 = [
1.0000     0.00189
0.9500     0.01210
0.9000     0.02172
0.8000     0.03935
0.7000     0.05496
0.6000     0.06845
0.5000     0.07941
0.4000     0.08705
0.3000     0.09003
0.2500     0.08912
0.2000     0.08606
0.1500     0.08018
0.1000     0.07024
0.0750     0.06300
0.0500     0.05332
0.0250     0.03922
0.0125     0.02841
0.0000     0.00000
0.0125     -0.02841
0.0250     -0.03922
0.0500     -0.05332
0.0750     -0.06300
0.1000     -0.07024
0.1500     -0.08018
0.2000     -0.08606
0.2500     -0.08912
0.3000     -0.09003
0.4000     -0.08705
0.5000     -0.07941
0.6000     -0.06845
0.7000     -0.05496
0.8000     -0.03935
0.9000     -0.02172
0.9500     -0.01210
1.0000     -0.00189
];
x_wing = NACA0018(:,1)' - 1;
y_wing = NACA0018(:,2)';
phi = -10/180*pi;
tmp = [cos(phi),-sin(phi);sin(phi),cos(phi)]*[x_wing; y_wing];
x_wing = tmp(1,:);
y_wing = tmp(2,:);

figure; pcolor(grid.X, grid.Y, reshape(Pref(30,:),size(X)))
shading interp; colormap(jet); axis equal; colorbar; caxis([-1 1]);
hold on; fill(x_wing, y_wing, 'k'); xlim([-1 2.4]); ylim([-0.5 0.5]);
title('REF');

figure; pcolor(grid.X, grid.Y, reshape(Prec(30,:),size(X)))
shading interp; colormap(jet); axis equal; colorbar; caxis([-1 1]);
hold on; fill(x_wing, y_wing, 'k'); xlim([-1 2.4]); ylim([-0.5 0.5]);
title('EPOD');

figure; pcolor(grid.X, grid.Y, reshape(Ptay(30,:),size(X)))
shading interp; colormap(jet); axis equal; colorbar; caxis([-1 1]);
hold on; fill(x_wing, y_wing, 'k'); xlim([-1 2.4]); ylim([-0.5 0.5]);
title('TH');

figure; pcolor(grid.X, grid.Y, reshape(Ppro(30,:),size(X)))
shading interp; colormap(jet); axis equal; colorbar; caxis([-1 1]);
hold on; fill(x_wing, y_wing, 'k'); xlim([-1 2.4]); ylim([-0.5 0.5]);
title('PRO');

%% functions
function V = Photo2world(V, Sample_rate, Mm_per_px_ratio)
% change velocity from photo to real world coordinate
V = V.*Sample_rate.*Mm_per_px_ratio.*1e-3;
end
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

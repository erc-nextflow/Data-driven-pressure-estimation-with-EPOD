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

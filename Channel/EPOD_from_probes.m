% this codes loads the fields and computes the POD
% Woii junwei.chen@uc3m.es
% v2 210517
% v3 210518
% v4 210518 use all probes
% v5 210521 fix bugs
% v6 210607
clear

Nsnap = 1:6400;
InputRootFields = './Fields/Field_';
InputRootProbes = './Probes/Probes_';
if ispc % change to Windows styled path
    InputRootFields(InputRootFields == '/') = '\';
    InputRootProbes(InputRootProbes == '/') = '\';
end
SavePath = 'EPOD_save.mat';
ProbeSelected = 1:10; % 1:20 whole; 1:10 trailing edge; 11:20 wall;
NoiseRatio = 0.00;

cont=0;
p_Wait = waitbar(0, 'Loading Fields...');
load(sprintf('%s%06d.mat',InputRootFields,Nsnap(1)),'u','xb','yb');
[X,Y] = meshgrid(xb,yb);
U = zeros(numel(u),numel(Nsnap));
V = U;
for i = Nsnap
    load(sprintf('%s%06d.mat',InputRootFields,i),'u','v');
    cont = cont+1; 
    U(:,cont) = u(:);
    V(:,cont) = v(:);
    if (mod(cont,100) == 0)
        tmp = (i-Nsnap(1))/(Nsnap(end)-Nsnap(1));
        waitbar(tmp*0.15, p_Wait, ['Loading Fields...', ...
            num2str(i), '/', num2str(numel(Nsnap))]);
    end
end

% adding noise
U = U + NoiseRatio*randn(size(U));
V = V + NoiseRatio*randn(size(V));


Um = mean(U,2);
Vm = mean(V,2);
U = U - Um;
V = V - Vm;

waitbar(0.15, p_Wait, 'Performing POD of fields...');
[psiF,sigmaF,phiF] = svd([U' V'],'econ');
waitbar(0.5, p_Wait, 'Performing POD of fields...');

figure
subplot(2,2,1)
yyaxis left
semilogx(diag(sigmaF).^2./sum(diag(sigmaF).^2)*100,'linewidth',1.5)
yyaxis right
semilogx(cumsum(diag(sigmaF).^2)./sum(diag(sigmaF).^2)*100,'linewidth',1.5)
subplot(2,2,2)
pcolor(X,Y,reshape(Um,size(X)))
shading interp
colormap jet(16)
title('mean flow - u component')
subplot(2,2,3)
pcolor(X,Y,reshape(phiF(1:numel(Um),1),size(X)))
shading interp
colormap jet(16)
title('First mode - u component')
subplot(2,2,4)
pcolor(X,Y,reshape(phiF(1+numel(Um):end,1),size(X)))
shading interp
colormap jet(16)
title('First mode - v component')

clear xb yb u v tmp i cont

waitbar(0.6, p_Wait, 'Loading Probe data...');
load(sprintf('%s%06d.mat',InputRootProbes,Nsnap(1)),'UP','VP');
UP = UP(:, ProbeSelected);
Upr = zeros(numel(UP), numel(Nsnap));
Vpr = Upr;
cont = 0;
for i = Nsnap
    load(sprintf('%s%06d.mat',InputRootProbes,i),'UP','VP');
    cont = cont + 1;
    Upr(:, cont) = reshape(UP(:, ProbeSelected), [size(Upr,1),1]);
    Vpr(:, cont) = reshape(VP(:, ProbeSelected), [size(Vpr,1),1]);
    if (mod(cont,10) == 0)
        tmp = (i-Nsnap(1))/(Nsnap(end)-Nsnap(1));
        waitbar(tmp*0.15+0.6, p_Wait, ['Loading Probe data...', ...
            num2str(i), '/', num2str(numel(Nsnap))]);
    end
end


% adding noise
Upr = Upr + NoiseRatio*randn(size(Upr));
Vpr = Vpr + NoiseRatio*randn(size(Vpr));

Um_pr = mean(Upr,2);
Vm_pr = mean(Vpr,2);
Upr = Upr - Um_pr;
Vpr = Vpr - Vm_pr;

waitbar(0.75, p_Wait, 'Performing POD of probes...');
[psiP,sigmaP,phiP] = svd([Upr' Vpr'],'econ');
waitbar(0.9, p_Wait, 'Finishing...');
% figure; plot(diag(sigmaP(1:12,1:12))./sum(diag(sigmaP)), '-o');

Xi = psiP'*psiF;

figure;
% imagesc(abs(Xi))
imagesc(abs(Xi(1:200, 1:200))); axis equal;
title('correlation matrix');

% filter to Xi
% Xi(abs(Xi)<(3/(sqrt(rank(Xi)))))=0;
% Xi(abs(Xi)<(3/(sqrt(length(Xi)))))=0;

MidMatrix = phiP/sigmaP*Xi*sigmaF;
clear phiP sigmaP sigmaF
save(SavePath, 'MidMatrix','phiF','X','Y','Um','Vm');
% save(SavePath, 'Um_pr','Vm_pr','Upr','Vpr','flag_ProbeGroup','-append');

close(p_Wait);

%% mean error map using training data

aExt = [Upr;Vpr]'*MidMatrix;
Urecon = aExt*phiF' + [Um;Vm]';
Uerr = sqrt(mean(((U+Um)' - Urecon(:,1:numel(X))).^2));
Verr = sqrt(mean(((V+Vm)' - Urecon(:,numel(X)+1:numel(X)*2)).^2));

figure; pcolor(X, Y, reshape(Uerr, size(X)));
shading interp; colormap(jet); axis equal; colorbar; % caxis([0 0.25]);
xlim([0 1]); ylim([0 1]);
title('MEAN EPOD ERROR of U (TRAINING DATA)');
figure; pcolor(X, Y, reshape(Verr, size(X)));
shading interp; colormap(jet); axis equal; colorbar; % caxis([0 0.25]);
xlim([0 1]); ylim([0 1]);
title('MEAN EPOD ERROR of V (TRAINING DATA)');

% figure; plot(aExt(200:400, 1), 'o-');
% title('time coefficient (fitst order)');

%% testing

Nsnap = 10001:11200;

cont = 0;
load(sprintf('%s%06d.mat',InputRootFields,Nsnap(1)),'u');
U = zeros(numel(u),numel(Nsnap));
V = U;
for i = Nsnap
    load(sprintf('%s%06d.mat',InputRootFields,i),'u','v');
    cont = cont+1; 
    U(:,cont) = u(:);
    V(:,cont) = v(:);
end

load(sprintf('%s%06d.mat',InputRootProbes,Nsnap(1)),'UP','VP');
UP = UP(:, ProbeSelected);
Upr = zeros(numel(UP), numel(Nsnap));
Vpr = Upr;
cont = 0;
for i = Nsnap
    load(sprintf('%s%06d.mat',InputRootProbes,i),'UP','VP');
    cont = cont + 1;
    Upr(:, cont) = reshape(UP(:, ProbeSelected), [size(Upr,1),1]);
    Vpr(:, cont) = reshape(VP(:, ProbeSelected), [size(Vpr,1),1]);
end

Upr = Upr - Um_pr;
Vpr = Vpr - Vm_pr;
aExt = [Upr;Vpr]'*MidMatrix;
Urecon = aExt*phiF' + [Um;Vm]';
Uerr = sqrt(mean((U'-Urecon(:,1:numel(X))).^2));
Verr = sqrt(mean((V'-Urecon(:,numel(X)+1:numel(X)*2)).^2));

figure; pcolor(X, Y, reshape(Uerr, size(X)));
shading interp; colormap(jet); axis equal; colorbar; % caxis([0 0.25]);
xlim([0 1]); ylim([0 1]);
title('MEAN EPOD ERROR of U (TESTING DATA)');
figure; pcolor(X, Y, reshape(Verr, size(X)));
shading interp; colormap(jet); axis equal; colorbar; % caxis([0 0.25]);
xlim([0 1]); ylim([0 1]);
title('MEAN EPOD ERROR of V (TESTING DATA)');

function [PIV,PR]=CalcPOD(PIV,PR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PIV.psi PIV.sigma PIV.phi]=svd([PIV.utr' PIV.vtr'],'econ');
[PR.psi PR.sigma PR.phi]=svd([PR.utr PR.vtr],'econ');

figure
subplot(1,2,1)
semilogx(diag(PIV.sigma.^2)./sum(diag(PIV.sigma.^2)),'o-')
subplot(1,2,2)
semilogx(cumsum(diag(PIV.sigma.^2))./sum(diag(PIV.sigma.^2)),'o-')


figure
subplot(1,2,1)
semilogx(diag(PR.sigma.^2)./sum(diag(PR.sigma.^2)),'o-')
subplot(1,2,2)
semilogx(cumsum(diag(PR.sigma.^2))./sum(diag(PR.sigma.^2)),'o-')
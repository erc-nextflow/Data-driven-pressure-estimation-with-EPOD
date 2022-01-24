function [U,V,PIV,PR,Hv,Wv,Nsnap]=ReadInput(Root,NImg,ROI,FlagFilt,TR,prx,pry,tau,Flagmean)


cont=0;
for i=NImg
    fprintf('%05d',i);
    cont=cont+1;
    INPUT=load(sprintf('%s%06d.mat',Root,i),'U','V');
    if i==NImg(1)
        U=zeros(ROI(2)-ROI(1)+1,ROI(4)-ROI(3)+1,numel(NImg));
        V=U;
    end
    U(:,:,cont)=INPUT.U(ROI(1):ROI(2),ROI(3):ROI(4));
    V(:,:,cont)=INPUT.V(ROI(1):ROI(2),ROI(3):ROI(4));
    fprintf('\b\b\b\b\b')
end
[Hv,Wv,Nsnap]=size(U);
if strcmp(FlagFilt,'YES')
    U = savitzkyGolay3D_rle_coupling(Hv,Wv,cont,U,7,7,7,3);
    V = savitzkyGolay3D_rle_coupling(Hv,Wv,cont,V,7,7,7,3);
end
%%%%%%%%%%%%%%%%%%%%%%
PIV.Utr=U(:,:,TR);
PIV.Vtr=V(:,:,TR);
PIV.Um=Flagmean*mean(PIV.Utr,3);
PIV.Vm=Flagmean*mean(PIV.Vtr,3);
PIV.utr=reshape(PIV.Utr-PIV.Um,[Hv*Wv numel(TR)]);
PIV.vtr=reshape(PIV.Vtr-PIV.Vm,[Hv*Wv numel(TR)]);


PR.Utr=zeros(numel(TR),tau*numel(pry));
PR.Vtr=PR.Utr;
cont=0;
for i=TR
    cont=cont+1;
    timevec=i:i+tau-1;
    PR.Utr(cont,:)=reshape(squeeze(U(pry,prx,timevec))',[1 numel(pry)*tau]);
    PR.Vtr(cont,:)=reshape(squeeze(V(pry,prx,timevec))',[1 numel(pry)*tau]);
end
PR.Um=mean(PR.Utr,1);
PR.Vm=mean(PR.Vtr,1);
PR.utr=PR.Utr-Flagmean*PR.Um;
PR.vtr=PR.Vtr-Flagmean*PR.Vm;


figure(1)
subplot(2,2,1)
imagesc(PIV.Um)
title('Mean crosswise velocity')
caxis([-1 1])
subplot(2,2,2)
imagesc(PIV.Vm)
caxis([12 18])
title('Mean streamwise velocity')
subplot(2,2,3)
imagesc(std(PIV.Utr,1,3))
caxis([0 5])
title('STD crosswise  velocity')
subplot(2,2,4)
imagesc(std(PIV.Vtr,1,3))
title('STD streamwise velocity')
caxis([0 5])
for i=1:4
    subplot(2,2,i)
    axis equal
    colormap jet(32)
    colorbar
end
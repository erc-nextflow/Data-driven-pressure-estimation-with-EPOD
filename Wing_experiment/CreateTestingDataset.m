function [PIV,PR]=CreateTestingDataset(U,V,TE,PIV,PR,Hv,Wv,tau,prx,pry,Flagmean)

PIV.ute=reshape(U(:,:,TE)-Flagmean*PIV.Um,[Hv*Wv numel(TE)]);
PIV.vte=reshape(V(:,:,TE)-Flagmean*PIV.Vm,[Hv*Wv numel(TE)]);

PR.Ute=zeros(numel(TE),tau*numel(pry));
PR.Vte=PR.Ute;
cont=0;
for i=TE
    cont=cont+1;
    timevec=i:i+tau-1;
    PR.Ute(cont,:)=reshape(squeeze(U(pry,prx,timevec))',[1 numel(pry)*tau]);
    PR.Vte(cont,:)=reshape(squeeze(V(pry,prx,timevec))',[1 numel(pry)*tau]);
end
PR.ute=PR.Ute-Flagmean*PR.Um;
PR.vte=PR.Vte-Flagmean*PR.Vm;

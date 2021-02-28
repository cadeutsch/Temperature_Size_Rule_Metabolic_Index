
%% Load climate data

load ../MAT/WOA
load ../MAT/cmip5_MM_dO2_dT.mat

z4d=permute(repmat(WOA.z',[360 1 180 12]),[1 3 2 4]);
o2fut=max(nanmean(WOA.o2,4)+dO2,0);
tmpfut=nanmean(WOA.temp,4)+dT;
po2fut=O2pressure(o2fut,tmpfut,tmpfut*0+35,z4d(:,:,:,1));
dpo2=po2fut-nanmean(WOA.po2,4);
po2a=nanmean(WOA.po2,4);

Iz=find(WOA.z==100);
plot(dT(:,:,Iz),dpo2(:,:,Iz),'k*')

%%

erange=-0.5:0.01:-0.05;
Erange=0:0.01:1;
[geps,gEo]=meshgrid(erange,Erange);

zmax=500;
dTz=dT(:,:,WOA.z<zmax);
dpo2z=dpo2(:,:,WOA.z<zmax);
po2z=po2a(:,:,WOA.z<zmax);
dTave=nanmean(dTz(:));
dpo2ave=nanmean(dpo2z(:))/nanmean(po2z(:));

dBBe=((1-f)./geps).*(gEo/(kb*(T0+Tref)^2)*dTave + dpo2ave);       
dBBe(dBBe<-1)=nan;        

%%

figure(1); clf reset

% lower panel
subplot('Position',[0.25 0.1 0.44 0.4]); 
[c,h]=contourf(erange,Erange,100*dBBe,-100:10:-10); clabel(c,h,'FontSize',14); %shading flat;caxis([0 1]);colorbar
xlabel('Allometric exponent (\epsilon)');ylabel('Temperature sensitivity (E_o)')

% upper panel
subplot('Position',[0.25 0.55 0.5 0.4]); 
Eo0=0.4; eps0=-0.3;
dBBe=(1-f)*(Eo0/(kb*(T0+Tref)^2)*dT - dpo2./po2a);
pcolor(100*nanmean(dBBe(:,:,1:10),3)'/eps0);shading flat;caxis([-30 30]); colorbar

[cmap]=cbrewer('div', 'RdBu', 10);
colormap(cmap)

print -djpeg TSE_Fig5.jpg
print(gcf,'-depsc2', '-painters', 'TSE_Fig5.eps');

%%


global ParInt

ParInt.B=B;
ParInt.epsilon=-epsilon;
ParInt.f=nanmedian(fSp);
ParInt.dpdT=0;

Tint=0:.1:10;
clear Bclim
Erange=[-0 0.1 0.4 0.7];
Brange=[1e-10 1e-4 1 1e2];
for jj=1:length(Erange)
    for ii=1:length(Brange)
    ParInt.Eo=Erange(jj);
    [t,y]=ode45(@dBdT_derivs,Tint,Brange(ii)); 
    Bclim(ii,jj,:)=(y/y(1)-1)*100;
    end
end


pcol='bcmrk';
ethresh=80;
figure(1); clf reset

xfill=[0 10 10 0 0];
yfill=[-100 -100 -80 -80 -100];
cfill=[1 1 1]*0.6;

subplot('Position',[0.15 0.15 0.3 0.3]); % lower left
for ii=1:length(Brange)
h(ii)=plot(Tint,squeeze(Bclim(ii,find(Erange==0.4),:)),pcol(ii)); hold on;
end
legstr=strcat(strtrim(cellstr(num2str(log10(Brange)',2))));
legstr={'10^{-10}','10^{-4}','10^{0}','10^{2}'};
legend(h,legstr,'Location','SouthWest','FontSize',14);text(0.4,-45,'Body mass [g]','FontSize',11)
%legend(legstr,'Location','East','FontSize',14);text(7.2,-30,'log_{10}(B) [g]','FontSize',12)
ylim([-100 10]);xlabel('Temperature Change (^oC)'); ylabel('Size Change (%)')
%plot(Tint,Tint*0-ethresh,'k--');text(1,-5-ethresh,'Extinction threshold?','FontSize',12)
fill(xfill,yfill,cfill,'FaceAlpha',0.4); text(6,-95,'Extinction?','FontSize',12)

subplot('Position',[0.55 0.15 0.3 0.3]); % lower right
for jj=1:length(Erange)
plot(Tint,squeeze(Bclim(find(Brange==100),jj,:)),pcol(jj)); hold on;
end
%legstr=strcat('E_o=',strtrim(cellstr(num2str(Erange',2))));
legstr=strcat(strtrim(cellstr(num2str(Erange',2))));
%legend(legstr,'Location','SouthWest','FontSize',14); text(0.5,-58,'E_o [eV]','FontSize',12)
legend(legstr,'Location','East','FontSize',14); text(7,-22,'Temp. Sens. [eV]','FontSize',10)
ylim([-100 10]);xlabel('Temperature Change (^oC)'); ylabel('Size Change (%)')
%plot(Tint,Tint*0-ethresh,'k--');text(1,-5-ethresh,'Extinction threshold?','FontSize',12)
fill(xfill,yfill,cfill,'FaceAlpha',0.4); text(0.2,-95,'Extinction?','FontSize',12)

% upper panel
subplot('Position',[0.15 0.55 0.73 0.4]); 
Eo0=0.4; eps0=-0.3;
dBBe=(1-f)*(Eo0/(kb*(T0+Tref)^2)*dT - dpo2./po2a);
pcolor(WOA.x,WOA.y,100*nanmean(dBBe(:,:,1:10),3)'/eps0);shading flat;caxis([-30 30]); 
h = colorbar;
set(get(h,'label'),'string','Size change (%)');
xlabel('Longitude');ylabel('Latitude')

[cmap]=cbrewer('div', 'RdBu', 10);
colormap(cmap)

print -djpeg TSE_Fig5.jpg
print(gcf,'-depsc2', '-painters', 'TSE_Fig5.eps');

%%

figure(2); clf reset
Eo0=0.4; eps0=-0.3;
subplot(211); 
dBBe=(1-f)*(Eo0/(kb*(T0+Tref)^2)*dT - 0*dpo2./po2a);
pcolor(100*nanmean(dBBe(:,:,1:10),3)'/eps0);shading flat;caxis([-30 30]); colorbar
subplot(212); 
dBBe=(1-f)*(Eo0/(kb*(T0+Tref)^2)*dT*0 - dpo2./po2a);
pcolor(100*nanmean(dBBe(:,:,1:10),3)'/eps0);shading flat;caxis([-30 30]); colorbar
[cmap]=cbrewer('div', 'RdBu', 10);
colormap(cmap)

print -djpeg TSE_FigS4.jpg
print(gcf,'-depsc2', '-painters', 'TSE_FigS4.eps');

set(groot,'defaultLineLineWidth', 2)
set(groot,'DefaultAxesFontSize', 20)
set(groot,'DefaultLegendFontSize', 14)


% TSE Figure 1

tmp_eps=-0.36;
tmp_Eo=nanmedian(Eo);
tmp_Ao=nanmedian(Ao);

brange=logspace(-1,1);
[X,Y]=meshgrid(brange,trange);
Tstar=(1/kb)*(1./(Y+T0) - 1/(Tref+T0));
Phi = Patm*tmp_Ao*(X.^tmp_eps).*exp(tmp_Eo*Tstar);


figure(1); clf reset;

subplot('Position',[0.2 0.55 0.4 0.4]); % [left bottom width height]
%subplot('Position',[0.1 0.3 0.4 0.4]); % [left bottom width height]
surf(Y,log10(X),Phi,'FaceAlpha',1);caxis([0 10]); c=colorbar; hold on; %shading flat; hold on;
contour3(Y,log10(X),Phi,[1 3 6],'k')
set(gca,'Ydir','reverse'); set(gca,'Xdir','reverse')
zlim([0 10]); view([-14 19])
ylabel('Body Mass','Rotation',-50); xlabel('Temperature'); 
zlabel('Metabolic Index \Phi'); 
%c.Label.String = 'Metabolic Index \Phi';

subplot('Position',[0.2 0.1 0.32 0.32]); % [left bottom width height]
[c,h]=contour(Y,log10(X),Phi,0:2:20,'LineWidth',3); hold on; % clabel(c,h,'FontSize',16)  %annotation('arrow',[0.8 0.9],[0.5 0.5]);
ylabel('Body Mass (log_{10})'); xlabel('Temperature'); axis([0 30 -1 1])

print -djpeg TSE_Fig1.jpg
print(gcf,'-depsc2', '-painters', 'TSE_Fig1.eps');


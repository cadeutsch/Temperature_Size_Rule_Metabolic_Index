%% MS Fig3 - Supply and Demand

figure(3); clf reset;

lfsize=12;

xl=[-11.5 3.5];
yl1=[-12 6];
yl2=[-3 3];

subplot(311);
plot(log10(Met.Mass),log10(Met.Rate),'r*'); hold on;
Btmp=Par.P.bstat(:,2);
I=find(isfinite(Btmp + alphaSb) & log10(Btmp)>-2);
plot(log10(Btmp(I)),log10(alphaSb(I)),'b*'); hold on; 
plot(log10(Supply_mass),log10(Supply_Ah),'bo'); hold on;

h1=legend('Demand (\alpha_D)','Supply (\alpha_S)','Gill (area)','Location','NorthWest');
set(h1,'FontSize',lfsize+2)

x=log10(Supply_mass); y=log10(Supply_Ah);
[p,pint,reg,rint,stats] = regress(y,[x ones(size(x))]);
sig_area=p(1);sig_area_err=abs(diff(pint(1,:)));
plot([min(x) max(x)],[min(x) max(x)]*p(1)+p(2),'--','Color',[1 1 1]*0.5)
plot(Bg,polyval(Met.p,Bg),'--','Color',[1 1 1]*0.5)
Ifit=find(log10(alphaSb)<30 & log10(Btmp)>-30); %find(log10(Par.M.bstat(:,2))<-2 | log10(Par.M.bstat(:,2))>3 | log10(alphaDb)>2.2);
x=log10(Btmp(Ifit));y=log10(alphaSb(Ifit));
[p,pint,reg,rint,stats] = regress(y,[x ones(size(x))]);
sig_inter=p(1); sig_inter_err=abs(diff(pint(1,:)));
plot([-1.7 max(x)],[-1.7 max(x)]*p(1)+p(2),'--','Color',[1 1 1]*0.5)
%xlabel('log10(B)');ylabel('log10(alphaSb)')
%text(0.1,0.9,['\sigma_{phy} = ' num2str(sig_inter(1),2)],'Units','Normalized','FontSize',20,'Color','b')
xlabel('log_{10}(Mass) [g]'); ylabel('log_{10}(Rate)'); %ylabel('log10(\alpha_D, \alpha_S)');
text(-4.3,-10.5,'TSE Size Range','FontSize',20,'Color',[0.5 0.5 0.5])
plot([b1 b1],yl1,'Color',[0.5 0.5 0.5]);plot([b2 b2],yl1,'Color',[0.5 0.5 0.5]);

xlim(xl); ylim(yl1); 

%%

subplot(312);


del=Par.M.est(:,2)+1; 
Idel=find(Par.flag_del & del>0.5 & del~=0.75);
del=del(Idel);
del_B=Par.M.bstat(Idel,2);
del_ont=[del_Clarke];%; del];
del_ont_B=[Mass_Clarke_med];%; del_B ];

I=find(isfinite(alphaSb) & Btmp>1);
H=errorbarxy(mean(log10(del_ont_B)),mean(del_ont),std(log10(del_ont_B)),std(del_ont), {'r', 'r', 'r'}); hold on
H.hMain.LineStyle='none';
errorbarxy(nanmean(log10(Btmp(I))),sig_inter,nanstd(log10(Btmp(I))),sig_inter_err,{'b', 'b', 'b'});
H.hMain.LineStyle='none';
errorbarxy(nanmedian(log10(Supply_mass)),sig_area,nanstd(log10(Supply_mass)),sig_area_err,{'b', 'b', 'b'});
H.hMain.LineStyle='none';

l1=plot(Bg(2:end),delta,'r-'); hold on
l2=plot(Bg(2:end),sigma,'b-')
l3=plot(mean(log10(del_ont_B)),mean(del_ont),'kp','MarkerFaceColor','r','MarkerSize',20); hold on;

%plot(log10(nanmean(Btmp(I))),sig_inter,'kp','MarkerFaceColor','b','MarkerSize',15)
l4=plot(nanmean(log10(Btmp(I))),sig_inter,'kp','MarkerFaceColor','b','MarkerSize',20);
l5=plot(nanmedian(log10(Supply_mass)),0.71,'kd','MarkerFaceColor','b','MarkerSize',10);

h2=legend([l1 l2 l3 l4 l5],'\delta_{phy}','\sigma (Model)','\delta_{ont}','\sigma_{phy}','Gill Area','Location','NorthWest')
set(h2,'FontSize',lfsize)

plot([b1 b1],[-.5 1.5],'Color',[0.5 0.5 0.5]);plot([b2 b2],[-.5 1.5],'Color',[0.5 0.5 0.5]);
plot([-12 4],[0 0],'k--')
xlim(xl); ylim([-0.2 1.5]); 
% plot(log10(nanmean(Par.B(Iem))),nanmedian(del(Iem))+1,'r^')
% plot(log10(Par.B(Iem)),del(Iem)+1,'r.')


xlabel('log_{10}(Mass) [g]'); ylabel('Allometric Exponent')

%%

subplot(313)

l0=plot(Bg(2:end),-epsilon,'g-'); hold on;
plot([b1 b1],[-1.5 .5],'Color',[0.5 0.5 0.5]);plot([b2 b2],[-1.5 .5],'Color',[0.5 0.5 0.5]);
text(-4.3,-1.3,'TSE Size Range','FontSize',20,'Color',[0.5 0.5 0.5])

% Inter-specific
%Btmp=Par.P.bstat(:,2);
I=find(log10(Btmp)>0 & Ao>0);
x=log10(Btmp(I));y=log10(Ao(I));
[p,pint,reg,rint,stats] = regress(y,[x ones(size(x))]);
eps_inter=p(1);eps_inter_err=abs(diff(pint(1,:)));
errorbarxy(nanmean(log10(Btmp(I))),eps_inter,nanstd(log10(Btmp(I))),eps_inter_err,{'k', 'k', 'k'});
l1=plot(nanmean(log10(Btmp(I))),eps_inter,'kp','MarkerFaceColor','g','MarkerSize',20)
H.hMain.LineStyle='none';

H=errorbarxy(Beps,Eps,lerrx, uerry*2,  {'k', 'k', 'k'});
H.hMain.LineStyle='none';
l2=plot(Beps,Eps,'ko','MarkerFaceColor','g','MarkerSize',9); hold on;

% [cmap]=cbrewer('seq','Greens',length(Eps)+2);
% msize=80;
% for ii=1:length(Eps)
%     H=errorbarxy(Beps(ii),Eps(ii),lerrx(ii), uerry(ii)*2,  {'k', 'k', 'k'});
%     scatter(Beps(ii),Eps(ii),msize,cmap(ii,:),'filled','MarkerEdgeColor',[0 .5 .5],'LineWidth',1.5); hold on;
% end

% data point from Calanus finmarchicus (see calanus_TSE.m)
ii=find(strcmp(TSE.Sp,'Calanus finmarchicus'));
plot(log10(nanmedian(TSE.B(ii,:)*1e-3)),-0.61,'ko','MarkerFaceColor','g','MarkerEdgeColor',[1 1 1]*.5,'MarkerSize',9)

% Plot Unicellular data
rmin=0.25;
Icil=strcmp(TSEuni2.taxa,'Ciliate') & TSEuni2.r2>rmin;
xcil=nanmedian(log10(w2d*TSEuni2.massave(Icil))); ycil=-nanmean(TSEuni2.eps(Icil));
excil=nanstd(log10(w2d*TSEuni2.massave(Icil))); eycil=nanstd(TSEuni2.eps(Icil));
errorbarxy(xcil,ycil,excil,eycil,{'go-', 'g', 'g'}); hold on;
l3=plot(xcil,ycil,'k^','MarkerFaceColor','g','MarkerEdgeColor',[1 1 1]*.5,'MarkerSize',9)
%Eps_add(end+1)=ycil;

Iame=strcmp(TSEuni2.taxa,'Amoeba') & TSEuni2.r2>rmin;
xame=nanmedian(log10(w2d*TSEuni2.massave(Iame))); yame=-nanmean(TSEuni2.eps(Iame));
exame=nanstd(log10(w2d*TSEuni2.massave(Iame))); eyame=nanstd(TSEuni2.eps(Iame));
errorbarxy(xame,yame,exame,eyame,{'go-', 'g', 'g'}); hold on;
plot(xame,yame,'k^','MarkerFaceColor','g','MarkerEdgeColor',[1 1 1]*.5,'MarkerSize',9)
%Eps_add(end+2)=yame;


plot([-12 4],[0 0],'k--')
xlabel('log_{10}(Mass) [g]'); ylabel('Allometric Exponent')
xlim(xl); ylim([-1.5 0.2]); 

legstr={'\epsilon (Model)','\epsilon_{phy}','\epsilon_{ont}','\epsilon_{uni}'};
h3=legend([l0 l1 l2 l3],legstr,'Location','SouthEast'); set(h3,'FontSize',lfsize)

Itmp=9:10;%strfind(sp_eps,'Scie')
hh=plot(Beps(Itmp),Eps(Itmp),'ko:','MarkerFaceColor','g','MarkerSize',13); 

print -djpeg TSE_Fig3.jpg
print(gcf,'-depsc2', '-painters', 'TSE_Fig3.eps');


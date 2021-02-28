%% Manuscript Fig2

figure(2); clf reset;
wysiwyg

set(gcf,'DefaultAxesFontSize', 14)

%% Histograms

clear frq_*

dV=1; Vedges=0:dV:30; 
dM=0.3; Medges=-2:dM:4; 
%dM=10; Medges=0:dM:600;
dE=0.2; Eedges=-1:dE:2; 
dD=2; Dedges=-5:dD:20;
dPhidT=100*Eo/(kb*(T0+Tref)^2);
eedges=-.65:.1:.6;
fedges=0:.2:1;

logBedge=-2:1:3;
[n,~,BbinEps]=histcounts(Beps,logBedge);
[n,~,BbinEo]=histcounts(log10(Par.B(Ieo)),logBedge);
[n,~,Bbinf]=histcounts(log10(Par.B),logBedge);
%log10(Par.B)

clear frq_*
%utaxP=phynames; %unique(Par.categ(Iem));
for i=1:length(logBedge)-1
    I=find(BbinEps==i);
    frq_Eps(:,i) = histc(Eps(I),eedges);
    I=find(BbinEo==i);
    frq_dPhidT(:,i) = histc(dPhidT(I),Dedges);
    I=find(Bbinf==i);
    frq_f(:,i) = histc(fSp(I),fedges);
    lgd{i}=[num2str(logBedge(i),2) ' to ' num2str(logBedge(i+1),2)];
%     frq_Eo(:,i) = histc(Eo(strcmp(utaxP{i},Par.categ))),Eedges);
%     frq_dPhidT(:,i) = histc((dPhidT(strcmp(utaxP{i},Par.categ))),Dedges);
end

[nanmean(-tse_sig(:)) nanstd(-tse_sig(:)) prctile(-tse_sig(:),[25 75])]
[nanmean(dPhidT) nanstd(dPhidT) prctile(dPhidT,[25 75])]

%% Plotting

[cmap]=cbrewer('seq','YlOrRd',size(frq_dPhidT,2));
colormap(cmap)
msize=100;
fsize=12;

figure(1); clf reset
subplot(321); hold off;
bar(Dedges+dD/2,frq_dPhidT,'stacked')
ylabel('# Species'); xlabel('E_{o}/k_B T^2 [%/^oC]','FontSize',14)
axis([-5 20 0 20])
legend(lgd,'Location','East','FontSize',10)
text(0.8,0.77,'log_{10}(B)','units','Normalized','FontSize',14)
colormap(cmap)

subplot(324)
X=Phi_max(Ieo); Y=Phi_min(Ieo); I=isfinite(X+Y);
[p,pint,~,~,stats] = regress(Y(I),[X(I) ones(size(X(I)))]);
%plot(Phi_max,Phi_min,'r*'); hold on;
I=BbinEo>0;
scatter(X(I),Y(I),X(I)*0+70,cmap(max(1,BbinEo(I)),:),'filled'); hold on;
plot(2*[0:10],[0:10],'k:');xlim([0 20])
xlabel('\Phi_{max}'); ylabel('\Phi_{crit}')
plot([0 20],[0 20],'k--')
plot([0 20],[0 10],'k:')
axis([0 20 0 20])

subplot(323)
bar(fedges+mean(diff(fedges))/2,frq_f,'stacked')
xlabel(['Activity buffer (f)'],'FontSize',16); ylabel('# Species')
xlim([0 1]);

subplot(325)
bar(eedges+mean(diff(eedges))/2,frq_Eps,'stacked')
xlabel(['Allometric exponent (' char(949) ')'],'FontSize',16); ylabel('# Species')
axis([-.6 .6 0 7])
legend(lgd,'Location','East','FontSize',10)
text(0.8,0.77,'log_{10}(B)','units','Normalized','FontSize',14)

print -djpeg TSE_Fig2.jpg
print(gcf,'-depsc2', '-painters', 'TSE_Fig2.eps');

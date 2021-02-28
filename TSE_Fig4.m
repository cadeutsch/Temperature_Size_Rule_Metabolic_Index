%% MS Fig4 - Predicted TSE

rlim=[0.25 1];
filter_tse

figure(4); clf reset;

set(0,'defaultLineLineWidth', 2)
set(0,'DefaultAxesFontSize', 16)
set(0,'DefaultLegendFontSize', 14)

Eedges=-2:.1:2;
Eaxes=[-.5 2 0 13]; 
Eaxes2=[-.5 2 0 45];

%%

tEo=(1-f)*100*Erange1/(kb*(T0+Tref)^2)/eps_pre;
tsr_edges=-20:0.2:20;
tse_fit=tEo(2:end)*nan;
for i=1:length(tEo)-1
    I=find(tEo>tsr_edges(i) & tEo<=tsr_edges(i+1));
    frq=Efit1(I);
    tse_fit(i)=sum(frq(:));
end
I=tse_fit>0;
tse_fit=tse_fit(I);
tsr_edges=tsr_edges(I);

%%

figure(4); clf reset

iAedges=[-30:2.5:5]; xl=[-20 5];
clear frq_tse
Hcateg='SZ';

if strcmp(Hcateg,'FM')
    [frq_tse(:,1),~]=histc(tse_sigM(:),iAedges);
    [frq_tse(:,2),~]=histc(tse_sigF(:),iAedges);
    disp('KS-test: Marine vs Freshwater:')
    [H,P]=kstest2(tse_sigM(:),tse_sigF(:))
    lgd={'Marine','Freshwater'};
    cmap2=cmap(2:4,:);
elseif strcmp(Hcateg,'SZ')
    logBedge=-6:1:0;
    [n,~,BbinTse]=histcounts(log10(tse_b(Isig)),logBedge);
    for i=1:length(logBedge)-1
        I=find(BbinTse==i);
        tmp=tse_sig(I,:);
        frq_tse(:,i) = histc(tmp(:),iAedges);
        lgd{i}=[num2str(logBedge(i),2) ' to ' num2str(logBedge(i+1),2)];
    end
    cmap2=cbrewer('seq','YlGnBu',length(logBedge)-1);
end

subplot(311); 
Idx=Isig; %find(tse_pval<1 & tse_r2>0.25);
tse_test=TSE.tse(Idx,:)*100;
%h1 = histogram(tse_test(:),iAedges); hold on;
%set(h1,'FaceColor',[1 1 1]*0.5,'EdgeColor','k','FaceAlpha',0.5);
h2 = bar(iAedges+mean(diff(iAedges))/2,frq_tse,'stacked'); hold on;
frq_tse_bar = histc(tse_bar(Isig),iAedges); 


plot(iAedges+mean(diff(iAedges))/2,frq_tse_bar*10,'r')
ylabel('# Samples, Species');xlabel('Observed TSE (%/^oC)')
colormap(cmap2);
xlim(xl)
hl=legend(lgd,'Location','East'); %legend(['All (p<1)',lgd])
set(hl,'FontSize',10)
text(0.88,0.8,'log_{10}(B)','units','Normalized','FontSize',12)

% bar(Vedges+dV/2,frq_Vh,'stacked')
% ylabel('# Species');xlabel('Hypoxia Vulnerability (V_h [atm])');
% axis([-0.05 0.2 0 21])
% hold on;
% h2a=legend(utaxP,'Position',[0.153 0.525 0.05 0.05] ); set(h2a,'FontSize',7)
% hAx(1)=gca;
% hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
% hold(hAx(2),'on')
% plot(hAx(2),Medges+dM/2,nansum(frq_alphaD,2),'k')
% axis([-2 3 0 42])
% xlabel('Metabolic Rate (log_{10}(\alpha_D))');
% h2b=legend('\alpha_D','Location','NorthEast'); set(h2b,'FontSize',10)

subplot(312); 
[npre,bin]=histc(tse_pred,iAedges); 
[nsig,bin]=histc(tse_sig(:),iAedges); 
[nall,bin]=histc(tse_test(:),iAedges); 
h1 = bar(iAedges+mean(diff(iAedges))/2,npre,'stacked'); hold on
sfac=max(npre)/max(nsig);
plot(iAedges+mean(diff(iAedges))/2,nsig*sfac,'k')
sfac=max(npre)/max(nall);
plot(iAedges+mean(diff(iAedges))/2,nall*sfac,'k--')
set(h1,'FaceColor','b','EdgeColor','k','FaceAlpha',0.5);
legend('Predicted','Observed')
ylabel('# Species');xlabel('TSE (%/^oC)')
xlim(xl);ylim([0 max(npre)+2])

subplot(313); 
%plot(log10(tse_b),tse_bar,'k*'); hold on
plot(log10(tse_b(Isig)),tse_bar(Isig),'ko'); hold on;

% p=polyfit(log10(tse_b(Isig)),tse_bar(Isig),1); 
% plot(x1,p(1)*x1+p(2),'k--')
%H=errorbar(log10(tse_b(Isig)),tse_bar(Isig),nanmin(tse_all(Isig,:),[],2),nanmax(tse_all(Isig,:),[],2),{'k', 'k', 'k'});
%plot(log10(tse_b(tse_pval<0.1)),tse_bar(tse_pval<0.1),'bo'); 
plot(log10b_tse,tse_slope,'b--'); 
legend('Data','Model','Location','NorthEast')
H=errorbar(log10(tse_b(Isig)),tse_bar(Isig),nanstd(tse_all(Isig,:),[],2)); %,{'k', 'k', 'k'});
H.LineStyle='none';
xlabel('log_{10}(B)');ylabel('TSE (%/degC)')
axis([-6 1 -15 5])

disp('KS-test: Obs vs Model:')
[H,P]=kstest2(tse_sig(:),tse_pred(:))

print -djpeg TSE_Fig4.jpg
print(gcf,'-depsc2', '-painters', 'TSE_Fig4.eps');



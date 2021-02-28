
clear all
close all

addpath(genpath('/Users/cdeutsch/Dropbox/code/m_files'))
addpath(genpath('/Users/cdeutsch/Dropbox/projects/MetHabClim/DATA'))
addpath(genpath('/Users/cdeutsch/Dropbox/projects/MetHabClim/code'))
addpath(genpath('/Users/cdeutsch/Dropbox/projects/MetHabClim/HeatWaves'))

set(groot,'defaultLineLineWidth', 2)
set(groot,'DefaultAxesFontSize', 20)
set(groot,'DefaultLegendFontSize', 14)
tfsize=16;
lfsize=16;

%% Set constants

global kb Rgas T0 Tref w2d rhoB Patm atm2kPa
global Glob

kb=8.6e-5; 
Rgas=8.3e-3;
T0=273.15;
Tref=15;
w2d=1; % wet to dry mass conversion
rhoB = 1; % density of biomass [g/cm3]
Patm=0.209;
atm2kPa=101.325;


%%  TSE data

load TSE

rlim=[0.25 1];
plim=0.05;
tse_pval=TSE.tse_stat(:,3);
tse_r2=TSE.tse_stat(:,1);
tse_b=nanmedian(TSE.B*1e-3,2); % grams
Isig=tse_r2>rlim(1) & tse_r2<rlim(2);
tse_sig=TSE.tse(Isig,:)*100;
tse_sigM=tse_sig(strcmp(TSE.dom(Isig),'Marine'),:);
tse_sigF=tse_sig(strcmp(TSE.dom(Isig),'Freshwater'),:);
tse_all=TSE.tse*100;
tse_bar=nanmean(tse_all,2);

%% Physiology Data

load('../MAT/Par.mat'); 
Eo=Par.P.est(:,3);
Ao=Par.P.est(:,1);
Iem=Par.Iem;
alphaD=Par.M.est(:,1); 
alphaS = alphaD.*Ao;
Btmp=Par.P.bstat(:,2);
alphaSb=alphaS.*Btmp;
alphaDb=alphaD.*Btmp;

Iem=Par.flag_em & isfinite(Par.M.est(:,1));
Rate_Brey=Par.M.est(Iem,1); % mg O2/day/mg body mass normalized to Tref=15
Mass_Brey=Par.M.bstat(Iem,2); % mg body mass
Rate_BreyW=Rate_Brey.*Mass_Brey; % mgO2/day
Par.Bm = Par.M.bstat(:,2); % median body mass for RMR
Par.Bp = Par.P.bstat(:,2); % median body mass for Pcrit
Par.B=Par.Bp;
Par.bminP=Par.P.bstat(:,3);
Par.bmaxP=Par.P.bstat(:,4);

% Eo 
Ieo=isfinite(Par.P.est(:,3));  Ieo(229:230)=0; % omit brittle-stars
Eo=Par.P.est(Ieo,3);
dE=0.5; 
E2=Eo+dE;
Erange2=[0:.01:2];
Epd=fitdist(E2,'lognormal');
Efit1 = pdf(Epd,Erange2);
Erange1=Erange2-dE;

[num,txt,raw]=xlsread('DeLong_BodyMass_vs_MetRate.xlsx'); % DeLong et al [2010]
Mass_DeLong=num(:,1); % gram
Rate_DeLong=num(:,2); % Watts
Rate_DeLong=Rate_DeLong*86400/14.1; % mgO2/day
Type_DeLong=txt(2:end,[1 4 5]); 
I=strcmp(txt(2:end,end),'Active');
Rate_DeLong(I)=Rate_DeLong(I)/2; % reduce active rates to resting (APPROX)

[num,txt,raw]=xlsread('Gillooly_BodyMass_vs_RBA.xlsx');
num(6,:)=nan; % obligate air-breather, omit
Supply_mass=num(:,1);
Supply_temp=num(:,2);
Supply_dz=num(:,3);
Supply_RBA=num(:,4);
Supply_Ah=Supply_RBA./Supply_dz;

[num,txt,raw]=xlsread('Clarke_BodyMass_vs_MetRate.xls'); % Clarke and Johnson [1999]
species2=txt(2:end,1);
% remove data from species already in Par
for i=1:length(species2)
    I=find(strcmp(Par.Species_unq,species2(i)));
    if ~isempty(I)
        num(i,:)=nan; 
    end
end
I=isfinite(num(:,1));
num=num(I,:);
species2=species2(I);

Temp_Clarke=num(:,2);
del_Clarke=num(:,4);
Rate_Clarke_O47g=num(:,7); % met rate @47g in O2 units
Rate_Clarke_med=num(:,10)*1e-3; % met rate @median size, in Watts
Mass_Clarke_med=num(:,8); % body mass of median size (g)
MetD2=num(:,7)./num(:,6)*1e3; % converted to ?mol O2/gram/hr

%% Biogeographic metrics for activity buffer (f)

load('../MAT/Pobis.mat')
phic=Pobis.phi_crit;
nmin=100;
Inum=Pobis.num_woa_unq<nmin;
Phi_max=Pobis.phi_wmax; Phi_max(Inum)=nan;
Phi_min=Pobis.phi_wmin; Phi_min(Inum)=nan;
fSp=1-Phi_min./Phi_max;

%% Epsilon estimates

load_Pcrit_vs_B

TSEuni2=load_TSR_unicell_2;


%% Allometry 

trange=0:30; 
B=logspace(-11,4,100);
Bg=log10(B); 
Bgcent=Bg(2:end)-mean(diff(Bg));

b1=min(log10(tse_b));
b2=max(log10(tse_b));

% Allometry - Demand
ord=3;mcut=-11;doff=0;
x=cat(1,log10(Mass_DeLong),[log10(Mass_Brey);  2.8]); % 4
y=cat(1,log10(Rate_DeLong),[log10(Rate_BreyW); 8+3 ]); % 3
p=polyfit(x(x>mcut),y(x>mcut),ord);
D_b = polyval(p,Bg);
delta=diff(D_b)./diff(Bg)+doff;

% Allometry - Supply
[n,t,r]=xlsread('lamellae.xlsx');
lam_d=n(:,2)*1e-4; % inter-lamellar distance (cm)
lam_l=n(:,3)*1e-4; % lamellar length (cm)
lam_h=lam_l; % lamellar height
lam_mass=n(:,1); % lam body mass
lam_area=lam_l.*lam_h; % lam area
I=log10(lam_mass)>4.5 | log10(lam_l)<-1.8;
lam_l(I)=nan; lam_h(I)=nan; %Lam.d(I)=nan;

% Model
% area of exchange elements
X=log10(lam_mass); Y=log10(lam_area); I=isfinite(X+Y);
[p,pint,~,~,stats] = regress(Y(I),[X(I) ones(size(X(I)))]);
a1=p(1); a0=10^p(2);
Ae=a0*B.^a1;
% radius of equivalent spherical exchange elements
r=(Ae/(4*pi)).^(1/2); % cm
L=(B/(4/3*pi*rhoB/1000)).^(1/3);
% distance
X=log10(lam_mass); Y=log10(lam_d); I=isfinite(X+Y);
[p,pint,~,~,stats] = regress(Y(I),[X(I) ones(size(X(I)))]);
d1=p(1); d0=10^p(2);
% cm/s velocity through gill
u0=1e-2;u1=0.8; 
U=u0*B.^u1;
% number of exchange elements
n0=1;n1=0.06; 
Bgill=1;
ne=n0*(B/Bgill).^n1;
%ne(ne<1)=1;


%% Optimize S model

C=[1 0.619 0.412]; % Sherwood coef's [Kiorboe et al. 2001]

optim_type='mech'; % functional ('func') or mechanistic ('mech')

optimize_sigma

% epsilon vs B
sigma=Sup.sigma;
epsilon=delta-sigma;

v2struct(Sup);

%% Compute TSE terms

% activity buffer 
f=nanmedian(fSp);  

% epsilon for TSE body sizes
x1=linspace(Bg(1),Bg(end));
I=x1>b1&x1<b2;
log10b_tse=x1(I);
eps_tse=epsilon(I);
eps_pre = -median(eps_tse); %-0.41; %

% prediction
term_Eo=100*(1-f)*Eo/(kb*(T0+Tref)^2)/eps_pre;
tse_pred=term_Eo;% -BxTbar-dEbar; %+term_dPhi;%
tse_slope=1*nanmean(term_Eo).*(mean(eps_tse))./eps_tse;

%% Plots 

TSE_Fig1

TSE_Fig2

TSE_Fig3

TSE_Fig4

TSE_Fig5


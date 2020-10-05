function [g0,g1,c,psi,ppi]=gammas(p,ss)

% fixed parameters
alphaa=0.33;
betta=0.99;
lambda=0.36;
epsilonp=6;
sigma=2;
thetap=0.75;
thetaw=0.75;
chiH=240;
chiL=24;

% estimated parameters
gammaH=p(1);
gammaL=p(2);
eta=p(3);
epsilonHw=p(4);
epsilonLw=p(5);
phip=p(6);
phiH=p(7);
phiL=p(8);
phiy=p(9);
rhoR=p(10);
rhonu=p(11);
rhoa=p(12);
rhoz=p(13);
rhoaH=p(14);
rhoaL=p(15);
rhop=p(16);
sigmanu=p(17);
sigmaa=p(18);
sigmaz=p(19);
sigmaaH=p(20);
sigmaaL=p(21);
sigmap=p(22);

% steady-state values
yss=ss(1);
cHss=ss(2);
cLss=ss(3);
nss=ss(4);
nHss=ss(5);
nLss=ss(6);
rwss=ss(7);
rwHss=ss(8);
rwLss=ss(9);

% Composite parameters
rho=-log(betta);
xiH=log(chiH);
xiL=log(chiL);
muHw=log(epsilonHw/(epsilonHw-1));
muLw=log(epsilonLw/(epsilonLw-1));
xpbar=log(epsilonp/(epsilonp-1));
etan=(eta-1)/eta;
etaw=1-eta;
Thetatp=(1-thetap)*(1-betta*thetap)*(1-alphaa)/(thetap*(1-alphaa+alphaa*epsilonp));
ThetaHw=(1-thetaw)*(1-betta*thetaw)/(thetaw*(1+gammaH*epsilonHw));
ThetaLw=(1-thetaw)*(1-betta*thetaw)/(thetaw*(1+gammaL*epsilonLw));
GammaHc=lambda*exp(cHss-yss);
GammaHn=lambda*(exp(nHss-nss)^etan);
GammaHw=lambda*(exp(rwHss-rwss)^etaw);

% equation indices
eq_nkpc=1;
eq_mc=2;
eq_Hwi=3;
eq_Hwmup=4;
eq_Lwi=5;
eq_Lwmup=6;
eq_euler=7;
eq_focn=8;
eq_rwH=9;
eq_rwL=10;
eq_ycl=11;
eq_nindex=12;
eq_windex=13;
eq_prodfn=14;
eq_Lbc=15;
eq_taylor=16;
eq_uH=17;
eq_uL=18;
eq_nu=19;
eq_a=20;
eq_z=21;
eq_aH=22;
eq_aL=23;
eq_xpt=24;
eq_Epip=25;
eq_EpiHw=26;
eq_EpiLw=27;
eq_EcH=28;
eq_yl=29;

% Variable indices
y=1;
cH=2;
cL=3;
n=4;
nH=5;
nL=6;
rw=7;
rwH=8;
rwL=9;
pip=10;
mupt=11;
piHw=12;
muHwt=13;
piLw=14;
muLwt=15;
i=16;
nuu=17;
a=18;
z=19;
aH=20;
aL=21;
xpt=22;
uH=23;
uL=24;
EcH=25;
Epip=26;
EpiHw=27;
EpiLw=28;
yl=29;

% Shock indices
omeganu=1;
omegaa=2;
omegaz=3;
omegaaH=4;
omegaaL=5;
omegap=6;

% Expectation errors
n_pip=1;
n_piHw=2;
n_piLw=3;
n_cH=4;

% Summary
neqs=29;
nshocks=6;
neta=4;
g0=zeros(neqs,neqs);
g1=zeros(neqs,neqs);
c=zeros(neqs,1);
psi=zeros(neqs,nshocks);
ppi=zeros(neqs,neta);

% 1. New Keynesian Phillips Curve
g0(eq_nkpc,pip)=1;
g0(eq_nkpc,Epip)=-betta;
g0(eq_nkpc,mupt)=Thetatp;
g0(eq_nkpc,xpt)=-Thetatp;

% 2. Marginal cost equation for firm
g0(eq_mc,mupt)=1;
g0(eq_mc,rw)=1;
g0(eq_mc,y)=alphaa/(1-alphaa);
g0(eq_mc,a)=-1/(1-alphaa);
c(eq_mc,1)=log(1-alphaa);

% 3. Wage inflation equation for HS households
g0(eq_Hwi,piHw)=1;
g0(eq_Hwi,EpiHw)=-betta;
g0(eq_Hwi,muHwt)=ThetaHw;
c(eq_Hwi,1)=ThetaHw*muHw;

% 4. Average wage mark-up for HS households
g0(eq_Hwmup,muHwt)=1;
g0(eq_Hwmup,rwH)=-1;
g0(eq_Hwmup,cH)=sigma;
g0(eq_Hwmup,nH)=gammaH;
c(eq_Hwmup,1)=xiH;

% 5. Wage inflation equation for LS households
g0(eq_Lwi,piLw)=1;
g0(eq_Lwi,EpiLw)=-betta;
g0(eq_Lwi,muLwt)=ThetaLw;
c(eq_Lwi,1)=ThetaLw*muLw;

% 6. Average wage mark-up for LS households
g0(eq_Lwmup,muLwt)=1;
g0(eq_Lwmup,rwL)=-1;
g0(eq_Lwmup,cL)=sigma;
g0(eq_Lwmup,nL)=gammaL;
c(eq_Lwmup,1)=xiL;

% 7. Euler equation for HS households
g0(eq_euler,cH)=1;
g0(eq_euler,EcH)=-1;
g0(eq_euler,i)=1/sigma;
g0(eq_euler,Epip)=-1/sigma;
g0(eq_euler,z)=-(1-rhoz)/sigma;
c(eq_euler,1)=rho/sigma;

% 8. FOC firm optimization problem
g0(eq_focn,nH)=1;
g0(eq_focn,nL)=-1;
g0(eq_focn,rwH)=eta;
g0(eq_focn,rwL)=-eta;
g0(eq_focn,aH)=-(eta-1);
g0(eq_focn,aL)=(eta-1);

% 9. HS real wage identity
g0(eq_rwH,rwH)=1;
g0(eq_rwH,piHw)=-1;
g0(eq_rwH,pip)=1;
g1(eq_rwH,rwH)=1;

% 10. LS real wage identity
g0(eq_rwL,rwL)=1;
g0(eq_rwL,piLw)=-1;
g0(eq_rwL,pip)=1;
g1(eq_rwL,rwL)=1;

% 11. Output market clearing condition
g0(eq_ycl,y)=1;
g0(eq_ycl,cH)=-GammaHc;
g0(eq_ycl,cL)=-(1-GammaHc);
c(eq_ycl,1)=yss-GammaHc*cHss-(1-GammaHc)*cLss;

% 12. Labor index
g0(eq_nindex,n)=1;
g0(eq_nindex,nH)=-GammaHn;
g0(eq_nindex,nL)=-(1-GammaHn);
g0(eq_nindex,aH)=-GammaHn;
g0(eq_nindex,aL)=-(1-GammaHn);
c(eq_nindex,1)=nss-GammaHn*nHss-(1-GammaHn)*nLss;

% 13. Wage index
g0(eq_windex,rw)=1;
g0(eq_windex,rwH)=-GammaHw;
g0(eq_windex,rwL)=-(1-GammaHw);
g0(eq_windex,aH)=GammaHw;
g0(eq_windex,aL)=(1-GammaHw);
c(eq_windex,1)=rwss-GammaHw*rwHss-(1-GammaHw)*rwLss;

% 14. Production function
g0(eq_prodfn,y)=1;
g0(eq_prodfn,n)=-(1-alphaa);
g0(eq_prodfn,a)=-1;

% 15. LS household budget constraint
g0(eq_Lbc,cL)=1;
g0(eq_Lbc,rwL)=-1;
g0(eq_Lbc,nL)=-1;

% 16. Taylor rule
g0(eq_taylor,i)=1;
g0(eq_taylor,pip)=-(1-rhoR)*phip;
g0(eq_taylor,piHw)=-(1-rhoR)*phiH;
g0(eq_taylor,piLw)=-(1-rhoR)*phiL;
g0(eq_taylor,y)=-(1-rhoR)*phiy;
g0(eq_taylor,nuu)=-1;
g1(eq_taylor,i)=rhoR;
c(eq_taylor,1)=(1-rhoR)*(rho-phiy*yss);

% 17. HS unemployment definition
g0(eq_uH,uH)=gammaH;
g0(eq_uH,muHwt)=-1;

% 18. LS unemployment definition
g0(eq_uL,uL)=gammaL;
g0(eq_uL,muLwt)=-1;

% SHOCKS
% 19. Monetary policy shock
g0(eq_nu,nuu)=1;
g1(eq_nu,nuu)=rhonu;
psi(eq_nu,omeganu)=1;

% 20. Technology shock
g0(eq_a,a)=1;
g1(eq_a,a)=rhoa;
psi(eq_a,omegaa)=1;

% 21. Preference (discount rate) shock
g0(eq_z,z)=1;
g1(eq_z,z)=rhoz;
psi(eq_z,omegaz)=1;

% 22. HS productivity shock
g0(eq_aH,aH)=1;
g1(eq_aH,aH)=rhoaH;
psi(eq_aH,omegaaH)=1;

% 23. LS productivity shock
g0(eq_aL,aL)=1;
g1(eq_aL,aL)=rhoaL;
psi(eq_aL,omegaaL)=1;

% 24. Price mark-up shock
g0(eq_xpt,xpt)=1;
g1(eq_xpt,xpt)=rhop;
c(eq_xpt,1)=(1-rhop)*xpbar;
psi(eq_xpt,omegap)=1;

% EXPECTATION ERRORS
% 24. Price inflation
g0(eq_Epip,pip)=1;
g1(eq_Epip,Epip)=1;
ppi(eq_Epip,n_pip)=1;

% 25. HS wage inflation
g0(eq_EpiHw,piHw)=1;
g1(eq_EpiHw,EpiHw)=1;
ppi(eq_EpiHw,n_piHw)=1;

% 26. LS wage inflation
g0(eq_EpiLw,piLw)=1;
g1(eq_EpiLw,EpiLw)=1;
ppi(eq_EpiLw,n_piLw)=1;

% 27. High-skill consumption
g0(eq_EcH,cH)=1;
g1(eq_EcH,EcH)=1;
ppi(eq_EcH,n_cH)=1;

% For measurement equation
% 28. lag value output, for output growth
g0(eq_yl,yl)=1;
g1(eq_yl,y)=1;

end
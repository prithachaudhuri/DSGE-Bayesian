function [A,B,H,R,Se,PHI]=sysmat(GG1,impact,p)

% errors
sigmanu=p(17);
sigmaa=p(18);
sigmaz=p(19);
sigmaaH=p(20);
sigmaaL=p(21);
sigmap=p(22);

% observables
eq_y=1;
eq_pip=2;
eq_ffr=3;
eq_uH=4;
eq_uL=5;
eq_sp=6;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transition equation
% s(t)=PHI*s(t-1)+R*e(t), e~N(0,Se)

nshocks=size(impact,2);
PHI=GG1;
R=impact;

Se=zeros(nshocks,nshocks);
Se(omeganu,omeganu)=sigmanu^2;
Se(omegaa,omegaa)=sigmaa^2;
Se(omegaz,omegaz)=sigmaz^2;
Se(omegaaH,omegaaH)=sigmaaH^2;
Se(omegaaL,omegaaL)=sigmaaL^2;
Se(omegap,omegap)=sigmap^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measurement equation
% y(t)=A+B*s(t)+u(t), u~N(0,H)

nobs=6;
nstates=size(PHI,2);

A=zeros(nobs,1);
% A(eq_y,1)=ss(1);
% A(eq_ffr,1)=-log(0.99);
% A(eq_ns,1)=ss(4);
% A(eq_nu,1)=ss(5);
% A(eq_uns,1)=ss(6);
% A(eq_unu,1)=ss(7);
% A(eq_sp,1)=ss(8);

B=zeros(nobs,nstates);
B(eq_y,y)=1;
B(eq_y,yl)=-1;
B(eq_pip,pip)=1;
B(eq_ffr,i)=1;
B(eq_uH,uH)=1;
B(eq_uL,uL)=1;
B(eq_sp,rwH)=1;
B(eq_sp,rwL)=-1;

H=zeros(nshocks,nshocks);


end
clear
clc
close all

% parameters
alphaa=0.33;
lambda=0.36; % CPS 2016
epsilonp=6;
sigma=2;
chiH=240;
chiL=24;
gammaH=7;
gammaL=2;
eta=2;
epsilonHw=4;
epsilonLw=7;

% composite parameters
mup=log(epsilonp/(epsilonp-1));
xiH=log(chiH);
xiL=log(chiL);
muHw=log(epsilonHw/(epsilonHw-1));
muLw=log(epsilonLw/(epsilonLw-1));
etan=(eta-1)/eta;
etaw=1-eta;

guess=[.5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 ... 
        .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 .5 ...
        .5 .5 .5 .5 .5 .5 .5 .5 .5]';
options=optimset('Display','iter','Tolfun',1.0e-8, 'MaxIter', 1.0e6, 'MaxFunEval', 1.0e6);
[g,fval,exitflag]=fsolve(@steady_state_solver, guess, options, sigma, xiH, xiL, alphaa, lambda, gammaH, gammaL, eta, ...
    muHw, muLw, mup, etan, etaw);

y=g(1);
cH=g(5);
cL=g(9);
n=g(13);
nH=g(17);
nL=g(21);
rw=g(25);
rwH=g(29);
rwL=g(33);

N=exp(n);
NH=exp(nH);
NL=exp(nL);
RWH=exp(rwH);
RWL=exp(rwL);
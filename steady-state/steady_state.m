clear
clc
close all

% parameters
alphaa=0.33;
epsilonp=6;
sigma=2;
chiH=1.757;
chiL=1.757;
gammaH=7;
gammaL=2;
epsilonHw=4;
epsilonLw=7;
lambda=0.25;
eta=1.2;

etan=(eta-1)/eta;
etaw=1-eta;

% steady-state values
Nss=0.23522;
Yss=Nss^(1-alphaa);
RWss=((epsilonp-1)/epsilonp)*(1-alphaa)*((1/Yss)^(1/(1-alphaa)));

guess=[0.5 0.5 0.5 0.5 0.5 0.5]';
options=optimset('Display','iter','Tolfun',1.0e-8, 'MaxIter', 1.0e10, 'MaxFunEval', 1.0e6);
[g,fval,exitflat]=fsolve(@solver,guess,options,Nss,Yss,RWss,sigma,chiH,chiL,gammaH,gammaL,epsilonHw,epsilonLw,lambda,eta);

% lb=0*ones(1,6);
% ub=0.99*ones(1,6);
% guess=0.5*ones(1,6);
% options=optimoptions(@lsqnonlin,'Display','Iter','FunctionTolerance',1e-25,'MaxFunctionEvaluations',200000,'MaxIterations',300000,'OptimalityTolerance',1e-25,'StepTolerance',1e-25);
% [ss,resnorm,residual,exitflag,output]=lsqnonlin(@(g) sssolver(g,Nss,Yss,RWss,sigma,chiH,chiL,gammaH,gammaL,epsilonHw,epsilonLw,lambda,etan,etaw,eta),guess,lb,ub,options);


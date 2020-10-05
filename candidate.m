%==========================================================================
% Skill heterogeneity in an estimated DSGE model
% Simplified model: 6 shocks
% Constructing candidate density for MH algorithm
% Pritha Chaudhuri, April 2018
%==========================================================================

clear 
clc
close all
delete *.asv
dbstop if error % stops code and shows where there is error
rng('default')

l=path;
path('mfiles',path);
path('gensys',path);
path('csminwell',path);
path('datasets',path);

tic

disp('                                                                  ');
disp('   Bayesian Estimation of DSGE model: candidate density           ');
disp('                                                                  ');


% parameters to estimate, starting values
param=[2.5 1.5 2.5 4 7 2 0.125 0.125 0.125 1.75 1 4 2.5 1 2.5 4 0.25 1 0.5 0.5 0.5 0.5];

[fh,x,gh,H,itct,fcount,retcodeh]=csminwel('objective',param',eye(22),[],10^(-5),500);

Theta=x';

Theta(10)=exp(Theta(10))/(1+exp(Theta(10)));
Theta(11)=exp(Theta(11))/(1+exp(Theta(11)));
Theta(12)=exp(Theta(12))/(1+exp(Theta(12)));
Theta(13)=exp(Theta(13))/(1+exp(Theta(13)));
Theta(14)=exp(Theta(14))/(1+exp(Theta(14)));
Theta(15)=exp(Theta(15))/(1+exp(Theta(15)));
Theta(16)=exp(Theta(16))/(1+exp(Theta(16)));

mode=Theta;
Sigma=hessian(@objective_uncons,mode');
Sigma=pinv(Sigma);

save datasets/mh_candidate mode Sigma

toc
%==========================================================================
% Skill heterogeneity in an estimated DSGE model
% Simplified model: 6 shocks
% MH algorithm
% Pritha Chaudhuri, April 2018
% Last edited: 07/10/2018 (posterior draws)
%==========================================================================

clear 
clc
close all
dbstop if error % stops code and shows where there is error
rng(0,'twister')

l=path;
path('mfiles',path);
path('gensys',path);
path('csminwell',path);
path('datasets',path);

tic

disp('                                                                  ');
disp('   Bayesian Estimation of DSGE model: MH algorithm           ');
disp('                                                                  ');

load mh_candidate
Sigma=nearestSPD(Sigma);

nsim=100000;
disp('                                                                  ');
disp('                                                                  ');

c=0.1; % scaling parameter for Hessian
c0=0.1;
% nburn=2;
nburn=int32(0.5*nsim)+2; % # burnin periods to discard
Thetasim=zeros(nsim,22); % 22 parameters to estimate

% Initialize Markov chain
% draw from Multivariate Normal Dist centered at mode with variance
% c0*Sigma. For check on the loop, draws for rho's must be between [0,1]
ctr=0;
while ctr==0
    Thetac=mvnrnd(mode',c0*Sigma);
    ctr=prod(Thetac(10:16)<=1); % loop ends if first draw satisfies bounds for rho's
end
Thetasim(1,:)=Thetac; % initialize the first simulation as the first draw

accept=0; % to calculate acceptance rate
obj=loglike(Thetasim(1,:))+prior(Thetasim(1,:)); % value of objective function at initialization
counter=0; % for # simualtions
logposterior=obj*ones(nsim,1);

for i=1:nsim
    
    % candidate draws for parameter values
    Thetac=mvnrnd(Thetasim(i,:),c*Sigma);
    checkbounds=prod(Thetac(10:16)<=1);
    
    % check if draws for rho's are in [0,1]
    if checkbounds==1
        priorc=prior(Thetac); % prior @ candidate draw
        llikec=loglike(Thetac); % log-likelihood value @ candidate draw
        objc=priorc+llikec; % objective function value @ candidate draw
        
        % check if objective function has infinite value, then next period
        % simulation value same as this period value
        if objc==-Inf
            Thetasim(i+1,:)=Thetasim(i,:);
            logposterior(i+1)=obj;
            
        else % if objective function value not negative infinity, then proceed with MH algorithm
            alpha=min(1,exp(objc-obj));
            u=rand(1);
            
            if u<=alpha % accept candidate draw
                Thetasim(i+1,:)=Thetac;
                accept=accept+1;
                obj=objc;
                logposterior(i+1)=objc;
            else % reject draw, next period same as this period
                Thetasim(i+1,:)=Thetasim(i,:);
                logposterior(i+1)=obj;
            end
            
        end % end objc==-Inf
        
    else % if checkbounds !=1 for rho's, next period draw same as this period
        Thetasim(i+1,:)=Thetasim(i,:);
        logposterior(i+1)=obj;
        
    end % end checkbounds==1
    
    % calculate acceptance rate and update counter for simulations
    acceptancerate=accept/i;
    counter=counter+1;
    
    % display simulation loop, acceptance rate and recursive averages for
    % parameter draws every 500 simulations
    if counter==500
        disp('                                                                   ');
        disp(['                       Draw number:', num2str(i)                  ]);
        disp('                                                                   ');
        disp(['                 Acceptance rate:', num2str(acceptancerate)       ]);
        disp('                                                                   ');
        disp('                  Recursive averages                               ');
        disp(num2str(mean(Thetasim(1:i,:))));
        disp('                                                                   ');
        counter=0;
    end
    
end

Thetasim=Thetasim(nburn:end,:);
logposterior=logposterior(nburn:end);

save datasets/mhdraws Thetasim logposterior


% Calculating posterior moments
post_mean=zeros(22,1);
post_median=zeros(22,1);
hpd_interval=zeros(22,2);

for i=1:22
    [post_mean(i,1),post_median(i,1),hpd_interval(i,:)]=post_moments(Thetasim(:,i));
end

toc
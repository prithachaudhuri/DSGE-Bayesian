function [post_mean,post_median,interval]=post_moments(xx)
xx=sort(xx);

post_mean=mean(xx);
post_median=median(xx);

ndraws=length(xx);
hpd_draws=round((1-0.95)*ndraws);

if hpd_draws>2
    kk=zeros(hpd_draws,1);
    jj=ndraws-hpd_draws;
    for ii=1:hpd_draws
        kk(ii)=xx(jj)-xx(ii);
        jj=jj+1;
    end
    [kmin,idx]=min(kk);
    interval=[xx(idx) xx(idx)+kmin];
else
    interval=NaN(1,2);
end
    
end
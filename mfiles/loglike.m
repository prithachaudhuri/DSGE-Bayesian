function [llhd,eu]=loglike(p)

% steady-state
ss=[-0.530633458221459; 0.036180583257200; -1.090975613064632; -0.791990236151431; -0.874415950113948; -0.747076299036864; -0.321442345461108; -0.280229488489212; -0.343899314027754];

% matrices for gensys
[g0,g1,c,psi,ppi]=gammas(p,ss);

% solving model
[GG1,C,impact,fmat,fwt,ywt,gev,eu,loose]=gensys(g0,g1,c,psi,ppi);

if eu(1)+eu(2)==2
    % state-space representation for Kalman Filter
    [A,B,H,R,Se,PHI]=sysmat(GG1,impact,p);
    
    load usdata7907; 
    
    % Kalman filter evaluation of likelihood function
    [loglike,mhat,shat,Phat]=kalman(A,B,H,R,Se,PHI,dy,pipobs,iobs,uHobs,uLobs,spobs);
    
    llhd=sum(loglike);
    
else
    
    llhd=-10^15;
    
end

end
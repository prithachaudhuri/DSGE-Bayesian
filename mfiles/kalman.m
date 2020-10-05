function [liki,measurepredi,statepredi,varstatepredi]=kalman(A,B,H,R,Se,PHI,dy,pipobs,iobs,uHobs,uLobs,spobs)
T=length(dy);
l=6;
data=zeros(T,l);
data(:,1)=dy;
data(:,2)=pipobs;
data(:,3)=iobs;
data(:,4)=uHobs;
data(:,5)=uLobs;
data(:,6)=spobs;

[n,n]=size(PHI);
s=zeros(T+1,n);
P=zeros(T+1,n,n);
s(1,:)=zeros(n,1)';

a=pinv(eye(n*n)-kron(PHI,PHI))*reshape(R*Se*R',n*n,1);
P(1,:,:)=reshape(a,n,n);

sprime=zeros(n,1);
Pprime=zeros(n,n);
errorprediction=ones(T,l);
Varerrorprediction=ones(T,l,l);
liki=ones(T,1);
measurepredi=ones(T,l);

for i=1:T
    
    % updating step
    sprime=PHI*s(i,:)';
    Pprime=PHI*squeeze(P(i,:,:))*PHI'+R*Se*R';
    
    % prediction step
    yprediction=A+B*sprime;
    
    v=data(i,:)'-yprediction;
    
    F=B*Pprime*B'+H;
    
    kgain=Pprime*B'*pinv(F);
    s(i+1,:)=(sprime+kgain*v)';
    P(i+1,:,:)=Pprime-kgain*B*Pprime;
    errorprediction(i,:)=v';
    Varerrorprediction(i,:,:)=F;
    liki(i)=real(-0.5*l*log(2*pi)-0.5*log(det(F))-0.5*v'*pinv(F)*v);
    measurepredi(i,:)=data(i,:)-v';
end

statepredi=s;
varstatepredi=P;

end

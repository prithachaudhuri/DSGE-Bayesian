function [p]=prior(Theta)
 % prior from uniform dist, params 10-16
 lb=[0,0,0,0,0,0,0];
 ub=[1,1,1,1,1,1,1];
 
 P10=1/(ub(1)-lb(1));
 P11=1/(ub(2)-lb(2));
 P12=1/(ub(3)-lb(3));
 P13=1/(ub(4)-lb(4));
 P14=1/(ub(5)-lb(5));
 P15=1/(ub(6)-lb(6));
 P16=1/(ub(7)-lb(7));
 
 
 % prior from normal dist
 P1=normpdf(Theta(1),2.5,0.5);
 P2=normpdf(Theta(2),1.5,0.25);
 P3=normpdf(Theta(3),1.5,0.25); % mean before 2.5
 
 
 % prior from gamma dist
 para1=[5,5,1,0.2,0.2,0.2];
 para2=[0.5,0.5,0.25,0.1,0.1,0.1];
 
 b=para2.^2./para1;
 a=para1./b;
 P4=gampdf(Theta(4),a(1),b(1));
 P5=gampdf(Theta(5),a(2),b(2));
 P6=gampdf(Theta(6),a(3),b(3));
 P7=gampdf(Theta(7),a(4),b(4));
 P8=gampdf(Theta(8),a(5),b(5));
 P9=gampdf(Theta(9),a(6),b(6));
 
 
 % prior from inverse gamma dist
 P17=exp(lnpdfig(Theta(17),0.5,4));
 P18=exp(lnpdfig(Theta(18),1,4));
 P19=exp(lnpdfig(Theta(19),0.5,4));
 P20=exp(lnpdfig(Theta(20),0.5,4));
 P21=exp(lnpdfig(Theta(21),0.5,4));
 P22=exp(lnpdfig(Theta(22),0.5,4));
 
 f=P1*P2*P3*P4*P5*P6*P7*P8*P9*P10*P11*P12*P13*P14*P15*P16*P17*P18*P19*P20*P21*P22;
 
 p=log(f);


end
function y = lnpdfig(x,a,b)
% LNPDFIG(X,A,B)
%	calculates log INVGAMMA(A,B) at X

% 03/03/2002
% Sungbae An
y = log(2) - gammaln(b/2) + (b/2)*log(b*a^2/2) - ( (b+1)/2 )*log(x^2) - b*a^2/(2*x^2);
end
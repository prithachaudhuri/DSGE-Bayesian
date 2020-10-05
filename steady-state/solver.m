function [f]=solver(g,Nss,Yss,RWss,sigma,chiH,chiL,gammaH,gammaL,epsilonHw,epsilonLw,lambda,eta)
f1 = g(5)-chiH*(g(1)^sigma)*(g(3)^gammaH)*((epsilonHw-1)/epsilonHw);
f2 = g(6)-chiL*(g(2)^sigma)*(g(6)^gammaL)*((epsilonLw-1)/epsilonLw);
% f3 = g(2)-g(6)*g(4);
f3 = g(3)/g(4)-((g(5)/g(6))^(-eta));
f4 = lambda*g(1)+(1-lambda)*g(2)-Yss;
f5 = lambda*(g(3)^((eta-1)/eta))+(1-lambda)*(g(4)^((eta-1)/eta))-(Nss^((eta-1)/eta));
f6 = lambda*(g(5)^(1-eta))+(1-lambda)*(g(6)^(1-eta))-(RWss^(1-eta));
f = [f1;f2;f3;f4;f5;f6];
end
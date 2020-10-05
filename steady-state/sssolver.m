function [ss]=sssolver(g,Nss,Yss,RWss,sigma,chiH,chiL,gammaH,gammaL,epsilonHw,epsilonLw,lambda,etan,etaw,eta)
ss=[g(5)-(epsilonHw/(epsilonHw-1))*chiH*(g(1)^sigma)*(g(3)^gammaH);
    g(6)-(epsilonLw/(epsilonLw-1))*chiL*(g(2)^sigma)*(g(4)^gammaL);
    (g(3)/g(4))-(g(5)/g(6))^(-eta);
    lambda*g(1)+(1-lambda)*g(2)-Yss;
    lambda*(g(3)^etan)+(1-lambda)*(g(4)^etan)-(Nss^etan);
    lambda*(g(5)^etaw)+(1-lambda)*(g(6)^etaw)-(RWss^etaw)];
end
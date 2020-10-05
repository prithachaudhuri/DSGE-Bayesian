function obj=objective(Theta)

Theta(10)=exp(Theta(10))/(1+exp(Theta(10)));
Theta(11)=exp(Theta(11))/(1+exp(Theta(11)));
Theta(12)=exp(Theta(12))/(1+exp(Theta(12)));
Theta(13)=exp(Theta(13))/(1+exp(Theta(13)));
Theta(14)=exp(Theta(14))/(1+exp(Theta(14)));
Theta(15)=exp(Theta(15))/(1+exp(Theta(15)));
Theta(16)=exp(Theta(16))/(1+exp(Theta(16)));


logprior=prior(Theta);

[llike,eu]=loglike(Theta);

% if eu(1)+eu(2)==-4
%     obj=-10^15;
% else
%     logprior=prior(Theta);
%     obj=-(llike+logprior);
% end

obj=-(llike+logprior);

end

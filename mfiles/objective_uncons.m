function [objective]=objective_uncons(Theta)
logprior=prior(Theta);

if logprior==-Inf
    objective=-10^9;
else 
    llike=loglike(Theta);
    objective=-(llike+logprior);
end

end
function J = exp_decay(k5,t,pErk)
% Likelihood function for the fit of the dephosphorylation rate.

J = 0;

b = -exp(-k5*max(t))/(1-exp(-k5*max(t)));
a = 1-b;

try
    simu = a*exp(-k5 * t) + b;
    for r = 1:size(pErk,1)
        J = J + nansum((simu-pErk(r,:)).^2);
    end
catch
    J = Inf;
end
end

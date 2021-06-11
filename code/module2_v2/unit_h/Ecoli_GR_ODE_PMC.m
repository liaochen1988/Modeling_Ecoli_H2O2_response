function [dxdt,mu] = Ecoli_GR_ODE_PMC(t,x,nutr,hp)

%   t: time
%   x: variables
%   nutr: concentration of growth-limiting nutrient
%   hp: host cell parameters

Met         = x(1);
Ribosome    = x(2);
ppGpp       = x(3);

%   Peptide elongation rate
k_pep  = hp.kelong*Met/(hp.Km_met+Met);
k_aap  = k_pep*Ribosome*hp.f_act;

%   Nutrient uptake and amino acid synthesis
MetE   = (hp.Pt*(1-hp.phiQ)-hp.Lr*Ribosome)*hp.alpha_metE/hp.Le;
if (MetE < 0)
    MetE = 0;
end

%   ppGpp synthesis rate
JsPpGpp = hp.ksPpGpp*hp.Km_met/(hp.Km_met+Met);

%   Ribosome synthesis rate
JsR     = k_aap*(1-hp.phiQ)*(hp.Ki_ppGpp./(hp.Ki_ppGpp+ppGpp))/hp.Lr;

%   Degradation rate of ppGpp
JdPpGpp = hp.kdPpGpp*ppGpp;

%   Specific growth rate
mu = k_aap/hp.Pt;

%   Differential equations
%   they are written in form of (synthesis rate)-(active
%   degradation)-(dilution)
dxdt(1) = nutr*MetE -k_aap*hp.fmet  -mu*Met;
dxdt(2) = JsR       -mu*Ribosome;
dxdt(3) = JsPpGpp   -JdPpGpp            -mu*ppGpp;
dxdt    = dxdt';

end


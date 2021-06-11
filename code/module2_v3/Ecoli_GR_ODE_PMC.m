function [dxdt,mu] = Ecoli_GR_ODE_PMC(t,x,nutr,H2O2,GSH,GSSG,Grx,hp)

%   t: time
%   x: variables
%   nutr: concentration of growth-limiting nutrient
%   hp: host cell paraMetE_SHrs

Met         = x(1);
Metox       = x(2);
Ribosome    = x(3);
MetE_SH     = x(4);
MetE_SSG    = x(5);
ppGpp       = x(6);

%   Peptide elongation rate
k_pep  = hp.kelong*Met/(hp.Km_met+Met);
k_aap  = k_pep*Ribosome*hp.f_act;

%   ppGpp synthesis rate
JsPpGpp = hp.ksPpGpp*hp.Km_met/(hp.Km_met+Met);

%   Ribosome synthesis rate
JsR     = k_aap*(1-hp.phiQ)*(hp.Ki_ppGpp./(hp.Ki_ppGpp+ppGpp))/hp.Lr;

%   MetE_SH synthesis rate
JsE     = k_aap*(1-hp.phiQ)*(ppGpp./(hp.Ki_ppGpp+ppGpp))*hp.alpha_metE/hp.Le;

%   Methionine oxidation rate
JoxMet  = hp.kox_met*H2O2*Met;

%   Methionine reduction rate
JredMet = hp.kred_met*Metox*hp.Trx_red*hp.Msr/(hp.Km_metox*Metox+hp.Km_trx*hp.Trx_red+Metox*hp.Trx_red);

%   MetE oxidation rate
JoxMetE = hp.kox_MetE*GSSG*MetE_SH;

%   MetE reduction rate
JredMetE = hp.kred_MetE*GSH*MetE_SSG;

%   MetE deglutathionylation by Grx
JredMetE_grx = hp.kred_MetE_grx*MetE_SSG*GSH^2*Grx/(MetE_SSG+hp.h_gssg*GSSG+hp.h_gsh*GSH^2);

%   Degradation rate of ppGpp
JdPpGpp = hp.kdPpGpp*ppGpp;

%   Specific growth rate
mu = k_aap/hp.Pt;

%   Differential equations
%   they are written in form of (synthesis rate)-(active
%   degradation)-(dilution)
dxdt(1) = nutr*MetE_SH -k_aap*hp.fmet  -mu*Met    -JoxMet + JredMet;
dxdt(2) = JoxMet    -JredMet        -mu*Metox;
dxdt(3) = JsR       -mu*Ribosome;
dxdt(4) = JsE       -mu*MetE_SH     -JoxMetE + JredMetE + JredMetE_grx;
dxdt(5) = JoxMetE   -JredMetE -JredMetE_grx - mu*MetE_SSG;
dxdt(6) = JsPpGpp   -JdPpGpp            -mu*ppGpp;
dxdt    = dxdt';

end


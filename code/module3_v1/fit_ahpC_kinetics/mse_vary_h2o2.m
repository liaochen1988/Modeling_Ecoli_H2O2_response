function [dy,Hout,Hin_f,frac_oxyr_ox_pred,Hin_J,Jahp_pred] = ...
    mse_vary_h2o2(p,h2o2_aslund1999_invivo,frac_oxyr_ox_aslund1999_invivo,h2o2_seaver2001,Jahp_seaver2001)

k_ahp = 10^p(1);
Km_ahp = 10^p(2);

%% Experiment 1: in vivo fraction of oxidized oxyR as a function of H2O2
% Aslund1999PNAS

% Use values from SeaverImlay2001JB for 
% P (permeability)
% A (area) 
% k_met (H2O2 production rate by metabolism)
% k_kat (H2O2 degradation rate of catalase)
P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20; % mol/s
k_kat = 1.05e-15; % mol/s
Km_kat = 3.9e-3; % M
Ki_katG_h2o2 = 19.55/1e6; % M, obtained from independent fitting
n_ahp = 1;

Hout = 10.^[-9:0.05:-4]; % M, extracellular H2O2
Hin_f = zeros(size(Hout)); % intracellular H2O2

options = optimset('Display','off','TolX',1e-12); % show iterations
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin_f(i+1);
    end
    
    [Hin_f(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_kat,Km_kat,Ki_katG_h2o2,k_ahp,Km_ahp,n_ahp);
    while Hin_f(i)<0
        x0 = 10.^(-9+rand*5);
        [Hin_f(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_kat,Km_kat,Ki_katG_h2o2,k_ahp,Km_ahp,n_ahp);
    end
    assert(exitflag>0);
end

% the following values also saved in par_est_oxyr_kinetics.mat
n_oxyr = 1.3614;
Km_oxyr = 41.3179/1e6; % M
kon_oxyr_h2o2 = 9.9863;
kon_oxyr_grx = 1.4910e7;
koff_oxyr_grx = 1.5218e15;

% Bionumber
% Grx1 ~200: Grx2 ~3,200: Grx3 ~4,400 (ng/mg)
% Grx1: 600*1e-9/1e-3/(1/1e-12)/9685/3.23e-15 = 19.2 uM
% Grx2: 2500*1e-9/1e-3/(1/1e-12)/24350/3.23e-15 = 31.79 uM
% Grx3: 4500*1e-9/1e-3/(1/1e-12)/9137/3.23e-15 = 152.48 uM

% We decide to use the in vitro parameters
GSSG = 1e-4; % M
GSH = 0.025; % M
Grx = 10*1e-6; % M
kswitch =  kon_oxyr_h2o2*Hin_f.^n_oxyr./(Km_oxyr^n_oxyr+Hin_f.^n_oxyr);
frac_oxyr_ox_pred = (kswitch+kon_oxyr_grx*Grx*GSSG^4)./(kswitch+kon_oxyr_grx*Grx*GSSG^4+koff_oxyr_grx*Grx*GSH^8);
frac_oxyr_ox_pred2 = pchip(Hout,frac_oxyr_ox_pred,h2o2_aslund1999_invivo);

%% Experiments: H2O2 decomposition rate as a function of H2O2 concentration
% SeaverImlay2001JB
% HPI-, k_kat = 0

Hin_J = zeros(size(Hout));
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin_J(i+1);
    end
    
    [Hin_J(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,0,Km_kat,Ki_katG_h2o2,k_ahp,Km_ahp,n_ahp);
    while Hin_J(i)<0
        x0 = 10.^(-9+rand*5);
        [Hin_J(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,0,Km_kat,Ki_katG_h2o2,k_ahp,Km_ahp,n_ahp);
    end
    assert(exitflag>0);
end
Jahp_pred = k_ahp.*Hin_J./(Hin_J+Km_ahp);
Jahp_pred2 = pchip(Hout,Jahp_pred,h2o2_seaver2001);

%% calculate difference
dy1 = (frac_oxyr_ox_pred2-frac_oxyr_ox_aslund1999_invivo)/max(frac_oxyr_ox_aslund1999_invivo);
dy2 = Jahp_pred2/max(Jahp_pred2)-Jahp_seaver2001/max(Jahp_seaver2001);
dy = [dy1,dy2]';

end


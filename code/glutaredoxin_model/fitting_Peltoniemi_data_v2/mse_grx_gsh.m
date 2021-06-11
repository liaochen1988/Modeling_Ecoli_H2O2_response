function dy = mse_grx_gsh(p,Grx1,initR1,Peptide2,initR2,Gsh3,initR3,Gssg4,initR4)

k1p = 10^p(1);
k2p = 10^p(2);
k2n = 10^p(3);

% Experiment 1: GSH = 1 mM, Peptide = 5 uM, GSSG = 0 uM
Grx = Grx1*1e-3; % uM
Gsh = 1*1e3; % uM
Peptide = 5; % uM
Gssg = 0; % uM
ypred1 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

% Experiment 2: GSH = 1 mM, Grx1 = 20 nM, GSSG = 0 uM
Grx = 20*1e-3; % uM
Gsh = 1*1e3; % uM
Peptide = Peptide2; % uM
Gssg = 0; % uM
ypred2 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

% Experiment 3: Peptide = 5 uM, Grx1 = 20 nM, GSSG = 0 uM
Grx = 20*1e-3; % uM
Gsh = Gsh3*1e3; % uM
Peptide = 5; % uM
Gssg = 0; % uM
ypred3 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

% Experiment 4: Peptide = 5 uM, Grx1 = 20 nM, GSH = 1 mM
Grx = 20*1e-3; % uM
Gsh = 1*1e3; % uM
Peptide = 5; % uM
Gssg = Gssg4; % uM
ypred4 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

dy = [ypred1-initR1,ypred2-initR2,ypred3-initR3,ypred4-initR4]';

end


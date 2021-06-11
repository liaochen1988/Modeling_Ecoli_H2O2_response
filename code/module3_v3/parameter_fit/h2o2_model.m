function dy = h2o2_model(t, x, constN, KatG_mutant, alpha_max_katG, alpha_max_ahpCF, alpha_max_grx, vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF)

H2O2_e = x(1);
H2O2_i = x(2);
OxyR_ox = x(3);
if (KatG_mutant)
    KatG = 0;
else
    KatG = x(4);
end
AhpCF = x(5);
Grx = x(6);
N = x(7);

% Parameters
P  = 1.6e-3;        % cm/s
A  = 1.41e-7;       % cm^2
PA = P * A / 1e3;   % L/s
k_met = 4.5e-20;    % mol/s
Ve = 0.02;          % L

kcat_katG       = 1.63e4;   % 1/s
Km_katG         = 3.9/1e3;  % M

kcat_ahpCF      = 52.4;     % 1/s
Km_ahpCF        = 1.5e-7;   % M

n_OxyR          = 1.36;
Km_OxyR         = 41.32e-6; % M
kon_OxyR_h2o2   = 9.99;     % 1/s
kon_OxyR_grx    = 1.49e7;   % 1/s/M^5
koff_OxyR_grx   = 1.52e15;  % 1/s/M^9
OxyR_tot        = 1.96e-6;  % M

GSSG = 0.00237;   % M
GSH  = 0.0166;    % M
mu   = 0.9/3600;  % 1/s, Glucose and glucose 6-phosphate as carbon sources in extra- and intracellular growth of enteroinvasive Escherichia coli and Salmonella enterica
Vc = 0.1882 * exp(1.6028 * mu * 3600) * 1e-15;

% OxyR oxidation/reduction
OxyR_red = OxyR_tot - OxyR_ox;
f = OxyR_ox / OxyR_tot;
v_OxyR_oxi_h2o2 = kon_OxyR_h2o2 * H2O2_i ^ n_OxyR / (Km_OxyR ^ n_OxyR + H2O2_i ^ n_OxyR) * OxyR_red;
v_OxyR_oxi_grx  = kon_OxyR_grx * Grx * GSSG ^ 4 * OxyR_red;
v_OxyR_red_grx  = koff_OxyR_grx * Grx * GSH ^ 8 * OxyR_ox;

% H2O2 degradation
v_katG_h2o2_deg  = kcat_katG * KatG * H2O2_i / (Km_katG + H2O2_i);
v_ahpCF_h2o2_deg = kcat_ahpCF * AhpCF * H2O2_i / (Km_ahpCF + H2O2_i);

% Catalase/Peroxidase deactivation rate
v_katG_deactiv  = vi_katG * KatG * H2O2_i / (H2O2_i + Ki_katG);
v_ahpCF_deactiv = vi_ahpCF * AhpCF * H2O2_i / (H2O2_i + Ki_ahpCF);

dy = zeros(7,1);
dy(1) = -(H2O2_e - H2O2_i) * PA / Ve * N;
dy(2) = (k_met + (H2O2_e-H2O2_i) * PA) / Vc - v_katG_h2o2_deg - v_ahpCF_h2o2_deg - v_katG_deactiv - v_ahpCF_deactiv;
dy(3) = v_OxyR_oxi_h2o2 + v_OxyR_oxi_grx - v_OxyR_red_grx;
dy(4) = alpha_max_katG * f - v_katG_deactiv - KatG  * mu;
dy(5) = alpha_max_ahpCF * f - v_ahpCF_deactiv - AhpCF * mu;
dy(6) = alpha_max_grx * f                   - Grx   * mu;
dy(7) = constN * (N * mu);


end
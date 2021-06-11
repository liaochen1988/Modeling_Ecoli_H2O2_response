function dx = Kinetic_model(t,x,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,called_by_ode15s)

dx = zeros(61,1);

y = real(abs(x));
Cell_N = y(1);
GLC_ex = y(2);
ACT_ex = y(3);
G6P = y(4);
FBP = y(5);
PG3 = y(6);
PEP = y(7);
PYR = y(8);
ACoA = y(9);
ICT = y(10);
AKG = y(11);
MAL = y(12);
OAA = y(13);
GLX = y(14);
PG6 = y(15);
Ribu5 = y(16);
Xyl5 = y(17);
Rib5 = y(18);
S7P = y(19);
E4P = y(20);
GLCin = y(21);
EIIA = y(22);
EIIA_P = y(23);
cAMP = y(24);
CrpcAMP = y(25);
CraFBP = y(26);
PdhRPYR = y(27);
SUC = y(28);
LAC = y(29);
AcAld = y(30);
ETH = y(31);
FOR = y(32);
QH2 = y(33);
Q = y(34);
NADH = y(35);
SUC_ex = y(36);

% new variables
ATP = y(37);
NADPH = y(38);

MET = y(39);
MetSO = y(40);
ppGpp = y(41);
PSH = y(42);
PSS = y(43);
MetE_SH = y(44);
MetE_SSG = y(45);
RIB = y(46);

H2O2_ex = y(47);
H2O2 = y(48);

OxyR_ox = y(49);
OxyR_red = y(50);

KatE = y(51);
if (katE_mutant)
    KatE = 0;
end
KatG = y(52);
if (katG_mutant)
    KatG = 0;
end
AhpCF = y(53);
Grx = y(54);
Gor = y(55);

GSH = y(56);
GSSG = y(57);
Trx_ox = y(58);
Trx_red = y(59);

NADP = y(60);
NAD = y(61);

Parameter

%   DO% == [DO2]/[DO2]_star*100
O2 = (a*0.214375)/100.0;

%% Total concentrations of TF and NADH
% Crp
Crp  = Crp_total - CrpcAMP;
% Cra
Cra  = Cra_total - CraFBP;
% PdhR
PdhR = PdhR_total - PdhRPYR;
% % NAD
% NAD = NADH_NAD_total - NADH;

%% Effect of TF on gene expression
% cAMPCrp
TF_CrpcAMP = CrpcAMP/(CrpcAMP+K_CrpcAMP);
% Cra
TF_Cra = Cra/(Cra+K_Cra);
% PdhR
TF_PdhR = PdhR/(PdhR+K_PdhR);
% IclR
TF_IclR = (CrpcAMP^n_IclR)/(K_IclR^n_IclR+CrpcAMP^n_IclR);
% Fnr
if fnr_mutant,
    TF_Fnr = 0.0;
else
    TF_Fnr = ((k_O2*O2)^n_FNR)/(K_A_FNR_O2^n_FNR+(k_O2*O2)^n_FNR);
end
% ArcA
if arcA_mutant,
    TF_ArcA = 0.0;
else
    TF_ArcA = (Q^n_ArcA)/(K_A_ArcA_Q^n_ArcA+Q^n_ArcA);
end


%% Table S2.  The kinetic equations for intracellular metabolic fluxes.
% PTS
% PTS1
v_PTS1 = PTS1_k1 * PEP * EIIA - PTS1_km1 * PYR * EIIA_P;
% PTS4
v_PTS4 = (1.0+TF_CrpcAMP) * PTS4_Vmax * EIIA_P * GLC_ex / ((PTS4_KEIIA + EIIA_P)*(PTS4_Kglc + GLC_ex));
% Non-PTS
% NPTS
v_NPTS = (1.0+TF_CrpcAMP) * NPTS_Vmax * GLC_ex / (NPTS_Ks + (1.0+(EIIA/NPTS_Ki))*GLC_ex);
% Glycolysis
% Glk
v_Glk = ( Glk_Vmax * GLCin ) / ( (1.0+(G6P/Glk_KI)) * (Glk_Ks+GLCin) );
% Pfk
v_Pfk = (1.0-TF_Cra) * Pfk_Vmax * (G6P/Pfk_Kg6p) * (1+G6P/Pfk_Kg6p)^(Pfk_n-1) / ((1+G6P/Pfk_Kg6p)^Pfk_n+Pfk_L*(1+PEP/Pfk_Kpep)^Pfk_n);
% Emp
v_Emp1 = (1.0-TF_Cra) * (Emp_Vmax_f*FBP/Emp_Kfbp-Emp_Vmax_r*PG3/Emp_Kpg3) / (1+FBP/Emp_Kfbp+PG3/Emp_Kpg3);
% Eno
v_Emp2 = (1.0-TF_Cra) * (Eno_Vmax_f*PG3/Eno_Kpg3-Eno_Vmax_r*PEP/Eno_Kpep) / (1+PG3/Eno_Kpg3+PEP/Eno_Kpep);
% Pyk
v_Pyk = (1.0-TF_Cra) * Pyk_Vmax * (PEP/Pyk_Kpep) * (1+PEP/Pyk_Kpep)^(Pyk_n-1) / ((1+PEP/Pyk_Kpep)^Pyk_n+Pyk_L/(1+FBP/Pyk_Kfbp)^Pyk_n);
% PDH
v_PDH = ((1.0-TF_PdhR)*(1.0-TF_ArcA)) * PDH_Vmax * PYR / PDH_Kpyr*(1+PYR/PDH_Kpyr)^(PDH_n-1)/((1+PYR/PDH_Kpyr)^PDH_n+PDH_L*(1+GLX/PDH_Kglx+PYR/PDH_KpyrI)^PDH_n);
% PP pathway
% Rpe
v_Ru5P = Ru5P_Vmax * (Ribu5-(Xyl5/Ru5P_Keq));
% Rpi
v_R5PI = R5PI_Vmax * (Ribu5-(Rib5/R5PI_Keq));
% TktA
v_TKa = TKa_Vmax * (Rib5*Xyl5-((S7P*PG3)/TKa_Keq));
% TktB
v_TKb = TKb_Vmax * (Xyl5*E4P-((G6P*PG3)/TKb_Keq));
% Tal
v_TA = TA_Vmax * (PG3*S7P-((E4P*G6P)/TA_Keq));
% FermentatiV_env fluxes
% LDH
v_LDH = Vmax_LDH * (NADH^n_NADH_LDH/(Km_NADH_LDH^n_NADH_LDH + NADH^n_NADH_LDH)) * (PYR^n_PYR_LDH/(Km_PYR_LDH^n_PYR_LDH + PYR^n_PYR_LDH));
% Pfl
v_Pfl = ((TF_Fnr+TF_ArcA)*Vmaxf_Pfl*PYR*CoA)/((PYR*CoA+Kpyr_Pfl*CoA+Kcoa_Pfl*PYR))...
    - ((TF_Fnr+TF_ArcA)*Vmaxr_Pfl*ACoA*FOR)/(ACoA*FOR+Kacoa_Pfl*FOR+Kfor_Pfl*ACoA);
% ALDH and ADH
v_ALDH = (Vmax_ALDH/(Kacoa_ALDH*Knadh_ALDH))*(ACoA*NADH-((CoA*NAD*AcAld)/Keq_ALDH))/((1.0+(NAD/Knad_ALDH)+(NADH/Knadh_ALDH))*(1.0+(ACoA/Kacoa_ALDH)+(CoA/Kcoa_ALDH)+(AcAld/Kacald_ALDH)+((AcAld*CoA)/(Kacald_ALDH*Kcoa_ALDH))));
v_ADH = ((Vmax_ADH/(Kacald_ADH*Knadh_ADH))*(AcAld*NADH-((ETH*NAD)/Keq_ADH)))/((1.0+(NAD/Knad_ADH)+(NADH/Knadh_ADH))*(1.0+(AcAld/Kacald_ADH)+(ETH/Keth_ADH)));
% Acetate formation/consumption
% PTACK
v_PTACK = PTACK_Vmax * (ACoA/PTACK_Kacoa) * (1+ACoA/PTACK_Kacoa)^(PTACK_n-1) / ((1+ACoA/PTACK_Kacoa)^PTACK_n+PTACK_L/(1+PYR/PTACK_Kpyr)^PTACK_n);
% Acs
v_Acs = (1.0-TF_ArcA)*(1.0-TF_Fnr)*(1.0+TF_CrpcAMP) * Acs_Vmax * ACT_ex / (ACT_ex + Acs_Kact);
% TCA cycle
% CS
v_CS = ((1.0+TF_CrpcAMP)*(1.0-TF_ArcA)) * CS_Vmax * OAA * ACoA / ((1+AKG/CS_Kakg)*CS_Koaaacoa*CS_Kacoa + CS_Kacoa*OAA + (1+AKG/CS_Kakg)*CS_Koaa*ACoA + OAA*ACoA);
% ICDH
v_ICDH = ((1.0+TF_Cra)*(1.0-TF_ArcA)) * ICDH_Vmax * (ICT/ICDH_Kict) * (1+ICT/ICDH_Kict)^(ICDH_n-1) / ((1+ICT/ICDH_Kict)^ICDH_n+ICDH_L*(1+PEP/ICDH_Kpep)^ICDH_n) * NADP / (NADP + KmNADP_ICDH);
% KGDH
v_KGMAL = ((1.0+TF_CrpcAMP)*(1.0-TF_ArcA)*(1.0-TF_Fnr)) * KGMAL_Vmax * AKG / (AKG + KGMAL_Kakg);
% SDH
v_SDH = ((1.0-TF_Fnr)*(1.0-TF_ArcA)*(1.0+TF_CrpcAMP))*(Vmaxf_SDH*SUC/Ksuc_SDH)/(1.0+(SUC/Ksuc_SDH)+(MAL/Kfum_SDH));
% Frd
v_FRD = (TF_Fnr)*(Vmaxr_SDH*MAL/Kfum_FRD)/(1.0+(SUC/Ksuc_SDH)+(MAL/Kfum_FRD));
% MDH
MDH_numerator = Vmaxf_MDH*MAL*NAD/(K_mal_MDH*K_nad_MDH) - Vmaxr_MDH*OAA*NADH/(K_oaa_MDH*K_nadh_MDH);
MDH_denominator = (1.0+MAL/K_mal_MDH+OAA/K_oaa_MDH)*(1.0+NAD/K_nad_MDH+NADH/K_nadh_MDH);
v_MDH = (1.0+TF_CrpcAMP)*(MDH_numerator/MDH_denominator);
% Anaplerotic reaction
% Ppc
v_Ppc = Ppc_Vmax * (PEP/Ppc_Kpep) * (1+PEP/Ppc_Kpep)^(Ppc_n-1) / ((1+PEP/Ppc_Kpep)^Ppc_n+Ppc_L/(1+FBP/Ppc_Kfbp)^Ppc_n);
% Glyoxylate pathways
% Icl
v_Icl = ((1.0+TF_Cra)*(1.0-TF_IclR)) * Icl_Vmax * (ICT/Icl_Kict) * (1+ICT/Icl_Kict)^(Icl_n-1) / ((1+ICT/Icl_Kict)^Icl_n+Icl_L*(1+PEP/Icl_Kpep+PG3/Icl_Kpg3+AKG/Icl_Kakg)^Icl_n);
% MS
v_MS = ((1.0+TF_Cra)*(1.0-TF_IclR)) * MS_Vmax * GLX * ACoA / (MS_Kglxacoa * MS_Kacoa + MS_Kacoa * GLX + MS_Kglx * ACoA + GLX * ACoA);
% Gluconeogenesis
% Fbp
v_Fdp = (1.0+TF_Cra) * Fdp_Vmax * (FBP/Fdp_Kfbp) * (1+FBP/Fdp_Kfbp)^(Fdp_n-1) / ((1+FBP/Fdp_Kfbp)^Fdp_n+Fdp_L/(1+PEP/Fdp_Kpep)^Fdp_n);
% Pps
v_Pps = (1.0+TF_Cra) * Pps_Vmax * (PYR/Pps_Kpyr) * (1+PYR/Pps_Kpyr)^(Pps_n-1) / ((1+PYR/Pps_Kpyr)^Pps_n+Pps_L*(1+PEP/Pps_Kpep)^Pps_n);
% Pck
v_Pck = (1.0+TF_Cra) * Pck_Vmax * OAA / (OAA+Pck_Koaa*(1+PEP/Pck_Kpep));
% Mez
v_Mez = Mez_Vmax * (MAL/Mez_Kmal) * (1+MAL/Mez_Kmal)^(Mez_n-1) / ((1+MAL/Mez_Kmal)^Mez_n+Mez_L*(1+ACoA/Mez_Kacoa+cAMP/Mez_Kcamp)^Mez_n) * NADP / (NADP + KmNADP_Mez);
% SUC transport
v_SUC_trans = SUC_trans_Vmax * SUC / (SUC_trans_Km_suc+SUC);

%% Kinetics for Cya, degradation of cAMP, and association/dissociation of TF
% Cya
v_Cya = Cya_Vmax * EIIA_P / (EIIA_P + Cya_KEIIA);
% cAMPdegr
v_cAMPdegr = cAMPdegr_Vmax * cAMP / (cAMP + cAMPdegr_KcAMP);
% cAMP-Crp
v_TF_cAMPCrp = TF_cAMPCrp_scale * ((Crp+CrpcAMP)*cAMP^TF_cAMPCrp_n/(cAMP^TF_cAMPCrp_n+TF_cAMPCrp_Kcamp^TF_cAMPCrp_n)-CrpcAMP);
% Cra-FBP
v_TF_CraFBP = TF_CraFBP_scale * ((Cra+CraFBP)*FBP^TF_CraFBP_n/(FBP^TF_CraFBP_n+TF_CraFBP_Kfbp^TF_CraFBP_n)-CraFBP);
% PdhR-PYR
v_TF_PdhRPYR = TF_PdhRPYR_scale * ((PdhR+PdhRPYR)*PYR^TF_PdhRPYR_n/(PYR^TF_PdhRPYR_n+TF_PdhRPYR_Kpyr^TF_PdhRPYR_n)-PdhRPYR);

%% Electron transport chain
% Nuo
v_Nuo = v_Nuo_max * (1.0-TF_ArcA) * (1.0-TF_Fnr) * (Q/(Km_Nuo_Q+Q)) * (NADH/(Km_Nuo_e2H2+NADH));
% Ndh
v_Ndh = v_Ndh_max * (1.0-TF_Fnr) * (Q/(Km_Ndh_Q+Q)) * (NADH/(Km_Ndh_e2H2+NADH));
% Cyo
v_Cyo = v_Cyo_max * (1.0-TF_ArcA) * (1.0-TF_Fnr) * (k_O2*O2/(Km_Cyo_O2+k_O2*O2)) * (QH2/(Km_Cyo_QH2+QH2));
% Cyd
v_Cyd = v_Cyd_max * (1.0+TF_ArcA) * (1.0-TF_Fnr) * (k_O2*O2/(Km_Cyd_O2+k_O2*O2)) * (QH2/(Km_Cyd_QH2+QH2));

%% Rate equations for fermentatiV_env variables
% Specific ATP production rate
OP_NADH = (4.0*v_Nuo + 4.0*v_Cyo + 2.0*v_Cyd)/3.0;
v_ATP = OP_NADH + v_Emp2 + v_Pyk + v_PTACK + v_KGMAL - v_Pfk - v_Pck - v_Acs - v_Glk;
v_NADPH_drain = kNADPH_drain * NADPH;

% Met synthesis
v_Bio_MET = k_ATP * v_ATP * MetE_SH * f_ATP_Met;
% Amino acid synthesis
v_Bio_AA = k_elong * MET / (Km_met + MET) * RIB * f_rib_act;
% MetE synthesis
v_Bio_MetE = v_Bio_AA * (1 - phi_Q) * (ppGpp / (Ki_ppGpp + ppGpp)) * alpha_MetE / l_MetE;
% Ribosome synthesis
v_Bio_RIB = v_Bio_AA * (1 - phi_Q) * (Ki_ppGpp / (Ki_ppGpp + ppGpp)) / l_Rib;
% ppGpp synthesis
v_Bio_ppGpp = k_ppGpp_syn * Km_met / (Km_met + MET);
% ppGpp degradation
v_Deg_ppGpp = k_ppGpp_deg * ppGpp;

% Specific growth rate
mu = v_Bio_AA / (PSH * 1e6);

% Cell volume
V_cell = 0.1882 * exp(1.6028 * mu * 3600) * 1e-15; % L

% The following reaction rates haV_env units g/L/s
% Glucose uptake rate
v_GLC_uptake = Cell_N * (v_PTS4 + v_NPTS) * uc_gDCW_L * V_cell / V_env * 1e-6 * Mw_GLC;
% Acetate excretion rate
v_ACT_excretion = Cell_N * v_PTACK * uc_gDCW_L * V_cell / V_env * 1e-6 * Mw_ACT;
% Acetate uptake rate
v_ACT_uptake = Cell_N * v_Acs * uc_gDCW_L * V_cell / V_env * 1e-6 * Mw_ACT;
% Lactate excretion rate
v_LAC_excretion = Cell_N * v_LDH * uc_gDCW_L * V_cell / V_env * 1e-6 * Mw_LAC ;
% Formate excretion rate
v_FOR_excretion = Cell_N * v_Pfl * uc_gDCW_L * V_cell / V_env * 1e-6 * Mw_FOR;
% Ethanol excretion rate
v_ETH_excretion = Cell_N * v_ADH * uc_gDCW_L * V_cell / V_env * 1e-6 * Mw_ETH;
% Succinate excretion rate
v_SUC_excretion = Cell_N * v_SUC_trans * uc_gDCW_L * V_cell / V_env * 1e-6 * Mw_SUC;

% The following reaction rates have unit M/s
% H2O2 transport (M/s)
v_H2O2_trans = (H2O2_ex - H2O2) * 69.85 * 5;
v_H2O2_trans_ex = Cell_N * v_H2O2_trans * V_cell / V_env;

%% G6PDH and PGDH
% G6PDH
v_G6PDH = (1.0-TF_Cra) * (mu*1.2202e+04) * G6PDH_Vmax * G6P * NADP / ((G6P+G6PDH_Kg6p)*(1.0+(NADPH/G6PDH_Knadph_g6pinh))*(G6PDH_Knadp*(1.0+(NADPH/G6PDH_Knadph_nadphinh))+NADP));
% PGDH
v_PGDH = (PGDH_Vmax * (mu*1.0044e+04) * PG6 * NADP) / ((PG6+PGDH_K6pg)*(NADP+PGDH_Knadp*(1.0+(NADPH/PGDH_Knadphinh))*(1.0+(ATP/PGDH_Katpinh))));

%% Synthesis
% Synthesis of QH2
v_syn_QH2 = mu * p_syn_QH2;
% Biosynthesis
% G6P
v_Bio_G6P = k_Bio_G6P * mu * G6P + b_Bio_G6P * G6P;
% GAP/�`/2PG
v_Bio_PG3 = k_Bio_PG3 * mu * PG3 + b_Bio_PG3 * PG3;
% PEP
v_Bio_PEP = k_Bio_PEP * mu * PEP + b_Bio_PEP * PEP;
% PYR
v_Bio_PYR = k_Bio_PYR * mu * PYR + b_Bio_PYR * PYR;
% AcCoA
v_Bio_ACoA = k_Bio_ACoA * mu * ACoA + b_Bio_ACoA * ACoA;
% R5P
v_Bio_R5P = k_Bio_R5P * mu * Rib5 + b_Bio_R5P * Rib5;
% E4P
v_Bio_E4P = k_Bio_E4P * mu * E4P + b_Bio_E4P * E4P;
% aKG
v_Bio_AKG = k_Bio_AKG * mu * AKG + b_Bio_AKG * AKG;
% OAA
v_Bio_OAA = k_Bio_OAA * mu * OAA + b_Bio_OAA * OAA;

%% OxyR dependent gene expression
frac_oxyr = OxyR_ox / (OxyR_ox + OxyR_red);

% KatG
v_Bio_katG = 5.50e-9 * frac_oxyr * 4;
% AhpCF
v_Bio_ahpCF = 3.48e-9 * frac_oxyr * 2;
% Grx
v_Bio_grx = 5.39e-9 * frac_oxyr;
% Gor
v_Bio_gor = 5.00e-9 * frac_oxyr;

%% H2O2 production and degradation
% Production in accompany with oxygen consumption rate
v_H2O2_metab_prod = (v_Cyo + v_Cyd) / 2 * uc_gDCW_L * 1e-6 * H2O2_leaky_ratio; % M/s
% Spontaneous decay
v_Deg_auto = 9.17e-6 * H2O2;
v_Deg_auto_ex = 9.17e-6 * H2O2_ex;
% KatG
v_Deg_h2o2_katG = 1.63e4 * KatG * H2O2 / (3.9e-3 + H2O2);
v_katG_deactiv = 2e-3 * KatG * H2O2 / (H2O2 + 1e-5);
% v_katG_deactiv = 20 * KatG * H2O2;
% KatE
v_Deg_h2o2_katE = 2.67e4 * KatE * H2O2 / (2.0e-2 + H2O2);
% AhpCF
v_Deg_h2o2_ahpCF = 52.4 * AhpCF * H2O2 / (1.50e-07 + H2O2) * NADH / (1.41e-6 + NADH);
v_ahpCF_deactiv = 200 * AhpCF * H2O2;

%% Oxidation/Reduction
% Methionine
v_Met_oxi = 6.0e-3 * H2O2 * MET;
v_Met_red = 3.7 * MetSO * Trx_red * MsrA / (1e-5 * MetSO + 1.9e-3 * Trx_red + MetSO * Trx_red);
% Thiol proteome
v_PSH_oxi = 1e4 * H2O2 * PSH;
v_PSH_red = 2.1e6 * PSS * Trx_red;
% MetE
v_MetE_oxi = 0.230 * MetE_SH * GSSG * 0;
v_MetE_red = 0.165 * MetE_SSG * GSH * 0;
v_MetE_red2 = 5.16e6 * MetE_SSG * GSH ^ 2  / (MetE_SSG + 1.451 * GSSG + 10.098 * GSH ^ 2) * Grx * 0;
% OxyR
v_OxyR_oxi = 9.99 * H2O2 ^ 1.36 / (4.13e-5 ^ 1.36 + H2O2 ^ 1.36) * OxyR_red;
v_OxyR_oxi2 = 1.49e7 * OxyR_red * GSSG ^ 4 * Grx;
v_OxyR_red2 = 1.52e15 * OxyR_ox * GSH ^ 8 * Grx;
% GSH
v_GSH_oxi = 1.5e9 * H2O2 * GSH ^ 2;
v_GSH_red = 1e2 * GSSG * Trx_red;
v_GSH_red2 = (267 * GSSG  + 6.55e5 * GSSG ^ 2) * (NADPH * uc_gDCW_L * 1e-6) * Gor / ...
             (2.22e-5 * GSSG + 9.7e-5 * (NADPH * uc_gDCW_L * 1e-6) + GSSG * (NADPH * uc_gDCW_L * 1e-6) + 0.022 * GSSG ^ 2 + 3.9e3 * GSSG ^ 2 * (NADPH * uc_gDCW_L * 1e-6)); % Gor-mediated reduction
% Thioredoxin
v_Trx_red = 41.25 * Trx_ox * (NADPH * uc_gDCW_L * 1e-6) * TrxR / (4.6e-6 * (NADPH * uc_gDCW_L * 1e-6) + 1.7e-6 * Trx_ox + (NADPH * uc_gDCW_L * 1e-6) * Trx_ox);

%% Mass balance equations

if called_by_ode15s,
    dx(1,1)  = Cell_N * mu;                                           % OD
    dx(2,1)  = - v_GLC_uptake;                                    % GLC_ex
    dx(3,1)  = v_ACT_excretion - v_ACT_uptake;                    % ACT_ex
    dx(29,1) = v_LAC_excretion;                                   % LAC
    dx(32,1) = v_FOR_excretion;                                   % FOR
    dx(31,1) = v_ETH_excretion;                                   % ETH
    % Glycolysis
    dx(21,1) = (v_NPTS) - (v_Glk) - mu * GLCin;                                                                   % GLCin
    dx(4,1)  = (v_PTS4 + v_Fdp + v_TKb + v_TA + v_Glk) - (v_Pfk + v_G6PDH + v_Bio_G6P) - mu * G6P;                % G6P/F6P
    dx(5,1)  = (v_Pfk) - (v_Emp1 + v_Fdp) - mu * FBP;                                                             % FBP
    dx(6,1)  = (2.0*v_Emp1 + v_TKa + v_TKb) - (v_Emp2 + v_TA + v_Bio_PG3) - mu * PG3;                             % GAP/DHAP
    dx(7,1)  = (v_Emp2 + v_Pck + v_Pps) - (v_Pyk + v_Ppc + v_PTS1 + v_Bio_PEP) - mu * PEP;                        % PEP
    dx(8,1)  = (v_Pyk + v_Mez + v_PTS1) - (v_PDH + v_Pps + v_Bio_PYR + v_LDH + v_Pfl) - mu * PYR;                 % PYR
    dx(9,1) = (v_PDH + v_Acs + v_Pfl) - (v_CS + v_PTACK + v_MS + v_Bio_ACoA + v_ALDH) - mu * ACoA;                % ACoA
    dx(30,1) = (v_ALDH) - (v_ADH) - mu * AcAld;                                                                   % AcAld
    % PP pathway
    dx(15,1) = (v_G6PDH) - (v_PGDH) - mu * PG6;                                                                   % PGL/6PG
    dx(16,1) = (v_PGDH) - (v_Ru5P + v_R5PI) - mu * Ribu5;                                                         % Ribu5
    dx(18,1) = (v_R5PI) - (v_TKa + v_Bio_R5P) - mu * Rib5;                                                        % Rib5
    dx(17,1) = (v_Ru5P) - (v_TKa + v_TKb) - mu * Xyl5;                                                            % Xyl5
    dx(19,1) = (v_TKa) - (v_TA) - mu * S7P;                                                                       % S7P
    dx(20,1) = (v_TA) - (v_TKb + v_Bio_E4P) - mu * E4P;                                                           % E4P
    % TCA cycle
    dx(10,1) = (v_CS) - (v_ICDH + v_Icl) - mu * ICT;                                                              % CIT/ICIT
    dx(11,1) = (v_ICDH) - (v_KGMAL + v_Bio_AKG) - mu * AKG;                                                       % AKG
    dx(28,1) = (v_KGMAL + v_Icl + v_FRD) - (v_SDH + v_SUC_trans) - mu * SUC;                                      % SUC
    dx(12,1) = (v_SDH + v_MS) - (v_MDH + v_Mez + v_FRD) - mu * MAL;                                               % MAL
    dx(13,1) = (v_MDH +v_Ppc) - (v_CS + v_Pck + v_Bio_OAA) - mu * OAA;                                            % OAA
    dx(14,1) = (v_Icl) - (v_MS) - mu * GLX;                                                                       % GOX
    % PTS protein, cAMP, and TF
    dx(22,1) = (v_PTS4) - (v_PTS1);                                                                               % EIIA
    dx(23,1) = (v_PTS1) - (v_PTS4);                                                                               % EIIA_P
    dx(24,1) = (v_Cya) - (v_cAMPdegr) - mu * cAMP;                                                                % cAMP
    dx(25,1) = (v_TF_cAMPCrp) - mu * CrpcAMP;                                                                     % CrpcAMP
    dx(26,1) = (v_TF_CraFBP) - mu * CraFBP;                                                                       % CraFBP
    dx(27,1) = (v_TF_PdhRPYR) - mu * PdhRPYR;                                                                     % PdhRPYR
    % Electron transport chain
    dx(33,1) = (v_Nuo + v_Ndh + v_SDH) - (v_Cyo + v_Cyd) + v_syn_QH2 - mu * QH2;                                  % QH2
    dx(34,1) = (v_Cyo + v_Cyd) - (v_Nuo + v_Ndh + v_SDH) - mu * Q;                                                % Q
    dx(35,1) = (v_Emp2 + v_PDH + v_KGMAL + v_MDH) - (v_LDH + v_ALDH + v_ADH + v_Nuo + v_Ndh) - v_Deg_h2o2_ahpCF;  % NADH
    dx(61,1) = - dx(35,1);                                                                                        % NAD
    % Extracellular SUC conc.
    dx(36,1) = v_SUC_excretion;                                                                                   % SUC_ex
    % ATP and NADPH
    dx(37,1) = v_ATP - (v_Bio_AA * 6 + v_Bio_RIB * 4566 * 10) / uc_gDCW_L - 0.82 * ATP - mu * ATP;                % ATP
    dx(38,1) = v_G6PDH + v_PGDH + v_ICDH + v_Mez - v_GSH_red2 - v_Trx_red - v_NADPH_drain;                        % NADPH
    dx(60,1) = - dx(38,1);                                                                                        % NADP
    % Methionine
    dx(39,1) = v_Bio_MET - v_Bio_AA * f_met  - v_Met_oxi + v_Met_red - mu * MET;                                  % MET
    dx(40,1) = v_Met_oxi - v_Met_red - mu * MetSO;                                                                % MetSO
    % ppGpp
    dx(41,1) = v_Bio_ppGpp - v_Deg_ppGpp - mu * ppGpp;                                                            % ppGpp
    % Thiol proteome
    dx(42,1) = 0; % - v_PSH_oxi + v_PSH_red;                                                                           % PSH
    dx(43,1) = 0; % - dx(42,1);                                                                                        % PSS
    % MetE
    dx(44,1) = v_Bio_MetE - v_MetE_oxi + v_MetE_red + v_MetE_red2 - mu * MetE_SH;                                 % MetE_SH
    dx(45,1) = v_MetE_oxi - v_MetE_red - v_MetE_red2 - mu * MetE_SSG;                                             % MetE_SSG
    % Ribosome   
    dx(46,1) = v_Bio_RIB - mu * RIB;                                                                              % RIB
    % H2O2
    dx(47,1) = - v_H2O2_trans_ex - v_Deg_auto_ex;                                                                 % H2O2_ex
    dx(48,1) = v_H2O2_metab_prod + v_H2O2_trans - v_Deg_auto ...
               - v_Deg_h2o2_katG - v_Deg_h2o2_katE - v_Deg_h2o2_ahpCF ...
               - n_oxyr * v_OxyR_oxi - 0 * v_PSH_oxi - v_Met_oxi - 0 * v_GSH_oxi;                                 % H2O2
    % OxyR
    dx(49,1) = v_OxyR_oxi + v_OxyR_oxi2 - v_OxyR_red2;                                                            % OxyR_ox
    dx(50,1) = - dx(49,1);                                                                                        % OxyR_red
    % OxyR-dependent gene upregulation
    dx(51,1) = 0;                                                                                                 % KatE
    dx(52,1) = v_Bio_katG - v_katG_deactiv - (mu + log(2)/Protein_half_life) * KatG;                                                                            % KatG
    dx(53,1) = v_Bio_ahpCF - v_ahpCF_deactiv - (mu + log(2)/Protein_half_life) * AhpCF;                                                                          % AhpCF
    dx(54,1) = v_Bio_grx - (mu + log(2)/Protein_half_life) * Grx;                                                                              % Grx
    dx(55,1) = v_Bio_gor - (mu + log(2)/Protein_half_life) * Gor;                                                                              % Gor
    % Glutathione
    dx(56,1) = - v_GSH_oxi * 2 + v_GSH_red * 2 + v_GSH_red2 * 2 ...
               + v_MetE_oxi - v_MetE_red - v_MetE_red2 ...
               + v_OxyR_oxi2 * 8 - v_OxyR_red2 * 8;                                                               % GSH                                  
    dx(57,1) = v_GSH_oxi - v_GSH_red - v_GSH_red2 ...
               - v_MetE_oxi + v_MetE_red + v_MetE_red2 ...
               - v_OxyR_oxi2 * 4 + v_OxyR_red2 * 4;                                                               % GSSG                                                                                  % GSSG
    % Thioredoxin
    dx(58,1) = - v_Trx_red + v_Met_red + v_PSH_red;                                                               % Trx_ox
    dx(59,1) = v_Trx_red - v_Met_red - v_PSH_red;                                                                 % Trx_red

    %   If t<0, run single-cell simulation by fixing external variables
    if t < 0.0
        dx(1,1)  = 0.0;                                   % N
        dx(2,1)  = 0.0;                                   % GLC_ex (external glucose)
        dx(3,1)  = 0.0;                                   % ACT_ex (external acetate)
        dx(29,1) = 0.0;                                   % LAC (lactate)
        dx(32,1) = 0.0;                                   % FOR (formate)
        dx(31,1) = 0.0;                                   % ETH (ethanol)
        dx(36,1) = 0.0;                                   % SUC_ex (external succinate)
        dx(47,1) = 0.0;                                   % H2O2_ex (external H2O2)                           
    end
    
else
    dx = [mu, v_H2O2_metab_prod];
    %   dx = [TF_Fnr TF_ArcA v_Cyo v_Cyd v_PDH v_Pfl v_LDH v_ADH v_FRD NAD];
end


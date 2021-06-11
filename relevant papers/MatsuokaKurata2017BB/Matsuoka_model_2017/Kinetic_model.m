function dx = Kinetic_model(t,x,arcA_mutant,fnr_mutant,a,called_by_ode15s)

dx = zeros(36,1);

OD = x(1);
GLC_ex = x(2);
ACT_ex = x(3);
G6P = x(4);
FBP = x(5);
PG3 = x(6);
PEP = x(7);
PYR = x(8); 
ACoA = x(9);
ICT = x(10);
AKG = x(11);
MAL = x(12);
OAA = x(13);
GLX = x(14);
PG6 = x(15);
Ribu5 = x(16);
Xyl5 = x(17);
Rib5 = x(18);
S7P = x(19);
E4P = x(20);
GLCin = x(21);
EIIA = x(22);
EIIA_P = x(23);
cAMP = x(24);
CrpcAMP = x(25);
CraFBP = x(26);
PdhRPYR = x(27);
SUC = x(28);
LAC = x(29);
AcAld = x(30);
ETH = x(31);
FOR = x(32);
QH2 = x(33);
Q = x(34);
NADH = x(35);
SUC_ex = x(36);
    
Parameter

O2 = (a*0.214375)/100.0;

%% Total concentrations of TF and NADH
% Crp
   Crp  = Crp_total - CrpcAMP;
% Cra
   Cra  = Cra_total - CraFBP;
% PdhR
   PdhR = PdhR_total - PdhRPYR;
% NAD
   NAD = NADH_NAD_total - NADH;
   
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
 % Fermentative fluxes
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
    v_ICDH = ((1.0+TF_Cra)*(1.0-TF_ArcA)) * ICDH_Vmax * (ICT/ICDH_Kict) * (1+ICT/ICDH_Kict)^(ICDH_n-1) / ((1+ICT/ICDH_Kict)^ICDH_n+ICDH_L*(1+PEP/ICDH_Kpep)^ICDH_n);   
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
    v_Mez = Mez_Vmax * (MAL/Mez_Kmal) * (1+MAL/Mez_Kmal)^(Mez_n-1) / ((1+MAL/Mez_Kmal)^Mez_n+Mez_L*(1+ACoA/Mez_Kacoa+cAMP/Mez_Kcamp)^Mez_n);
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

%% Rate equations for fermentative variables
% Specific ATP production rate
   OP_NADH = (4.0*v_Nuo + 4.0*v_Cyo + 2.0*v_Cyd)/3.0;
   v_ATP = OP_NADH + v_Emp2 + v_Pyk + v_PTACK + v_KGMAL - v_Pfk - v_Pck - v_Acs - v_Glk;
   
% Specific growth rate
   mu =  k_ATP * v_ATP;

% Glucose uptake rate
   v_GLC_uptake = Mw_GLC * uc * OD * (v_PTS4 + v_NPTS);
% Acetate excretion rate
   v_ACT_excretion = Mw_ACT * uc * OD * v_PTACK;
% Acetate uptake rate
   v_ACT_uptake = Mw_ACT * uc * OD * v_Acs;
% Lactate excretion rate
   v_LAC_excretion = Mw_LAC * uc * OD * v_LDH ;
% Formate excretion rate
   v_FOR_excretion = Mw_FOR * uc * OD * v_Pfl;
% Ethanol excretion rate
   v_ETH_excretion = Mw_ETH * uc * OD * v_ADH;
% Succinate excretion rate
   v_SUC_excretion = Mw_SUC * uc * OD * v_SUC_trans;

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
  % GAP/Å`/2PG
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


%% Mass balance equations

if called_by_ode15s,
   dx(1,1)  = OD * mu;                                           % OD
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
   dx(35,1) = (v_Emp2 + v_PDH + v_KGMAL + v_MDH) - (v_LDH + v_ALDH + v_ADH + v_Nuo + v_Ndh) - mu * NADH;         % NADH   
 % Extracellular SUC conc.
   dx(36,1) = v_SUC_excretion;                                                                                   % SUC_ex
   
   if t < 0.0
   dx(1,1)  = 0.0;                                   % OD
   dx(2,1)  = 0.0;                                   % GLC_ex
   dx(3,1)  = 0.0;                                   % ACT_ex
   dx(29,1) = 0.0;                                   % LAC
   dx(32,1) = 0.0;                                   % FOR
   dx(31,1) = 0.0;                                   % ETH
   dx(36,1) = 0.0;                                   % SUC_ex
   end
   
else
   dx = [TF_Fnr TF_ArcA v_Cyo v_Cyd v_PDH v_Pfl v_LDH v_ADH v_FRD NAD];
end


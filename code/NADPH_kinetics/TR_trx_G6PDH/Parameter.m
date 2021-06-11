%% Parameters for fermentative variables
k_ATP = 0.00001;         % [gDCW/É mol]
Mw_GLC = 180.156000;
Mw_XYL = 150.130000;
Mw_ACT = 60.050000;
Mw_LAC = 90.08;
Mw_FOR = 46.025;
Mw_ETH = 46.07;
Mw_SUC = 118.09;
uc = 9.5E-07;            % [gDCW/É (OD).l]

%% Kinetic parameters for intracellular fluxes
% PTS
 % PTS1
   PTS1_k1 = 116.000000; % [gDCW/gProt.s]
   PTS1_km1 = 46.3;      % [gDCW/gProt.s]
 % PTS4
   PTS4_Vmax = 4.3766;   % [É mol/gDCW.s]
   PTS4_KEIIA = 0.0085;  % [gProt/gDCW]
   PTS4_Kglc = 0.001200; % [g/l]
% Non-PTS
 % NPTS
   NPTS_Vmax = 0.1670;   % [É mol/gDCW.s]
   NPTS_Ks = 0.01;       % [É mol/gDCW]
   NPTS_Ki = 5.0e-06;    % [É mol/gDCW]
% Glycolysis
 % Glk
   Glk_Vmax = 5.556;     % [É mol/gDCW.s]
   Glk_KI = 0.100000;    % [É mol/gDCW]
   Glk_Ks = 16.000000;   % [É mol/gDCW]
 % Pfk
   Pfk_Vmax = 92.510592; % [É mol/gDCW.s]
   Pfk_Kg6p = 0.022000;  % [É mol/gDCW]
   Pfk_Kpep = 0.138000;  % [É mol/gDCW]
   Pfk_L = 9.5e+07;      %
   Pfk_n = 4.000000;     %
 % Emp
   Emp_Vmax_f = 10.194683; % [É mol/gDCW.s]
   Emp_Vmax_r = 10.194683; % [É mol/gDCW.s]
   Emp_Kfbp = 5.920000;    % [É mol/gDCW]
   Emp_Kpg3 = 16.600000;   % [É mol/gDCW]
 % Eno
   Eno_Vmax_f = 10.194683; % [É mol/gDCW.s]
   Eno_Vmax_r = 10.194683; % [É mol/gDCW.s]
   Eno_Kpg3 = 4.760000;    % [É mol/gDCW]
   Eno_Kpep = 1.110000;    % [É mol/gDCW]
 % Pyk
   Pyk_Vmax = 49.501994;   % [É mol/gDCW.s]
   Pyk_Kpep = 5.000000;    % [É mol/gDCW]
   Pyk_Kfbp = 0.413000;    % [É mol/gDCW]
   Pyk_L = 100000.000000;  %
   Pyk_n = 4.000000;       %
 % PDH
   PDH_Vmax = 5.6667;      % [É mol/gDCW.s]
   PDH_Kpyr = 0.128000;    % [É mol/gDCW]
   PDH_Kglx = 0.218000;    % [É mol/gDCW]
   PDH_KpyrI = 0.231000;   % [É mol/gDCW]
   PDH_L = 3.400000;       %
   PDH_n = 2.650000;       %
% PP pathway
 % G6PDH
   G6PDH_Vmax = 2.4472*10000;              % [É mol/gDCW.s]
   G6PDH_Kg6p = 25.5319;             % [É mol/gDCW]
   G6PDH_Knadp = 0.0436*100;             % [É mol/gDCW]
   G6PDH_Knadph_g6pinh = 11.4007;    % [É mol/gDCW]
   G6PDH_Knadph_nadphinh = 0.0177/100;   % [É mol/gDCW]
 % PGDH
   PGDH_Vmax = 28.7808;              % [É mol/gDCW.s]
   PGDH_K6pg = 66.4894;              % [É mol/gDCW]
   PGDH_Knadp = 0.0897;              % [É mol/gDCW]
   PGDH_Knadphinh = 0.0245 * 10000;          % [É mol/gDCW]
   PGDH_Katpinh = 368.7943;          % [É mol/gDCW]
 % Rpe
   Ru5P_Vmax = 11.9486;              % [É mol/gDCW.s]
   Ru5P_Keq = 1.400000;              %
 % Rpi
   R5PI_Vmax = 8.5787;               % [É mol/gDCW.s]
   R5PI_Keq = 4.000000;              %
 % TktA
   TKa_Vmax = 16.7968;               % [É mol/gDCW.s]
   TKa_Keq = 1.200000;               %
 % TktB
   TKb_Vmax = 153.4726;              % [É mol/gDCW.s]
   TKb_Keq = 10.000000;              %
 % Tal
   TA_Vmax = 19.2760;                % [É mol/gDCW.s]
   TA_Keq = 1.050000;                %
% Fermentative fluxes
 % LDH
   Vmax_LDH = 51.7;                  % [É mol/gDCW.s]
   Km_NADH_LDH = 0.9;                % [É mol/gDCW]
   Km_PYR_LDH = 26.7;                % [É mol/gDCW]
   n_NADH_LDH = 3.0;                 %
   n_PYR_LDH = 2.0;                  %
 % Pfl
   Vmaxf_Pfl = 14.7;                 % [É mol/gDCW.s]
   Vmaxr_Pfl = 3.4;                  % [É mol/gDCW.s] 
   Kpyr_Pfl = 4.5;                   % [É mol/gDCW]
   Kcoa_Pfl = 0.7816;                % [É mol/gDCW]
   Kacoa_Pfl = 5.8621;               % [É mol/gDCW]
   Kfor_Pfl = 1.1276;                % [g/l]
 % ALDH
   Vmax_ALDH = 63.9;                 % [É mol/gDCW.s]
   Kacoa_ALDH = 0.0124;              % [É mol/gDCW]
   Knadh_ALDH = 0.443;               % [É mol/gDCW]  
   Knad_ALDH = 0.1418;               % [É mol/gDCW]
   Kcoa_ALDH = 0.0142;               % [É mol/gDCW]
   Kacald_ALDH = 17.7305;            % [É mol/gDCW]
   Keq_ALDH = 2.7730e+01;            % [É mol/gDCW]   
 % ADH
   Vmax_ADH = 57.97;                 % [É mol/gDCW.s]
   Kacald_ADH = 0.0532;              % [É mol/gDCW]
   Knadh_ADH = 2.1;                  % [É mol/gDCW]
   Knad_ADH = 0.1418;                % [É mol/gDCW]
   Keth_ADH = 7.1065e-04;            % [g/l]
   Keq_ADH = 2.1906e+04;             % [É mol/gDCW]
% Acetate formation/consumption
 % PTACK
   PTACK_Vmax = 2.2;                 % [É mol/gDCW.s]
   PTACK_Kacoa = 0.022000;           % [É mol/gDCW]
   PTACK_Kpyr = 0.022000;            % [É mol/gDCW]
   PTACK_L = 6.39e+05;               %
   PTACK_n = 2.000000;               %
 % Acs
   Acs_Vmax = 0.7294;                % [É mol/gDCW.s]
   Acs_Kact = 0.001;                 % [g/l]
% TCA cycle
 % CS
   CS_Vmax = 2.1154;                 % [É mol/gDCW.s]
   CS_Kacoa = 0.212000;              % [É mol/gDCW]
   CS_Koaa = 0.029000;               % [É mol/gDCW]
   CS_Koaaacoa = 0.029000;           % [É mol/gDCW]
   CS_Kakg = 0.630000;               % [É mol/gDCW]
 % ICDH
   ICDH_Vmax = 7.8873;               % [É mol/gDCW.s]
   ICDH_Kict = 0.000160;             % [É mol/gDCW]
   ICDH_Kpep = 0.334000;             % [É mol/gDCW]
   ICDH_L = 127.000000;              %
   ICDH_n = 2.000000;                %
 % KGDH
   KGMAL_Vmax = 1.1832;              % [É mol/gDCW.s]
   KGMAL_Kakg = 0.548000;            % [É mol/gDCW]
 % SDH
   Vmaxf_SDH = 1.8101;               % [É mol/gDCW.s]
   Ksuc_SDH = 0.2299;                % [É mol/gDCW]
   Kfum_SDH = 0.5747;                % [É mol/gDCW]
 % Frd
   Vmaxr_SDH = 2.5;                  % [É mol/gDCW.s]
   Kfum_FRD = 0.5747;                % [É mol/gDCW]
 % MDH
   Vmaxf_MDH = 12.6770;              % [É mol/gDCW.s]
   Vmaxr_MDH = 54.3318;              % [É mol/gDCW.s]
   K_mal_MDH = 29.8;                 % [É mol/gDCW]
   K_nad_MDH = 29.8851;              % [É mol/gDCW]
   K_oaa_MDH = 0.56322;              % [É mol/gDCW]
   K_nadh_MDH = 7.0115;              % [É mol/gDCW]
% Anaplerotic reaction
 % Ppc
   Ppc_Vmax = 1.0;                   % [É mol/gDCW.s]
   Ppc_Kpep = 0.048000;              % [É mol/gDCW]
   Ppc_Kfbp = 0.408000;              % [É mol/gDCW]
   Ppc_L = 5200000.000000;           %
   Ppc_n = 3.000000;                 %
% Glyoxylate pathways
 % Icl
   Icl_Vmax = 1.2324;                % [É mol/gDCW.s]
   Icl_Kict = 0.0022;                % [É mol/gDCW]
   Icl_Kpep = 0.055000;              % [É mol/gDCW]
   Icl_Kpg3 = 0.720000;              % [É mol/gDCW]
   Icl_Kakg = 0.827000;              % [É mol/gDCW]
   Icl_L = 5.01;                     %
   Icl_n = 4.000000;                 %
 % MS
   MS_Vmax = 2.1444;                 % [É mol/gDCW.s]
   MS_Kglx = 0.950000;               % [É mol/gDCW]
   MS_Kacoa = 0.755000;              % [É mol/gDCW]
   MS_Kglxacoa = 0.719000;           % [É mol/gDCW]
% Gluconeogenesis
 % Fbp
   Fdp_Vmax = 0.1140;                % [É mol/gDCW.s]
   Fdp_Kfbp = 0.003000;              % [É mol/gDCW]
   Fdp_Kpep = 0.300000;              % [É mol/gDCW]
   Fdp_L = 4000000.000000;           %
   Fdp_n = 4.000000;                 %
 % Pps
   Pps_Vmax = 0.0041;                % [É mol/gDCW.s]
   Pps_Kpyr = 0.001770;              % [É mol/gDCW]
   Pps_Kpep = 0.001000;              % [É mol/gDCW]
   Pps_L = 1.0e-79;                  %
   Pps_n = 2.000000;                 %
 % Pck
   Pck_Vmax = 0.0872;                % [É mol/gDCW.s]
   Pck_Koaa = 0.184000;              % [É mol/gDCW]
   Pck_Kpep = 1000.000000;           % [É mol/gDCW]
 % Mez
   Mez_Vmax = 0.707450;              % [É mol/gDCW.s]
   Mez_Kmal = 0.006240;              % [É mol/gDCW]
   Mez_Kacoa = 3.640000;             % [É mol/gDCW]
   Mez_Kcamp = 6.540000;             % [É mol/gDCW]
   Mez_L = 104000.000000;            %
   Mez_n = 1.330000;                 %
% SUC transport
   SUC_trans_Vmax = 0.5;             % [É mol/gDCW.s]
   SUC_trans_Km_suc = 16.0;          % [É mol/gDCW] 
% Biosynthetic pathway
   % G6P
     k_Bio_G6P = 554.4;              %
     b_Bio_G6P = 0.002;
   % GAP/Å`/2PG
     k_Bio_PG3 = 352.8000;           %
     b_Bio_PG3 = 0.002;
   % PEP
     k_Bio_PEP = 3.0456e+03;         %
     b_Bio_PEP = 0.002;
   % PYR
     k_Bio_PYR = 39.816;             %
     b_Bio_PYR = 0.002;
   % AcCoA
     k_Bio_ACoA = 135.36;            %
     b_Bio_ACoA = 0.035;
   % R5P
     k_Bio_R5P = 72.0;               %
     b_Bio_R5P = 0.002;
   % E4P
     k_Bio_E4P = 72.0;               %
     b_Bio_E4P = 0.002;
   % aKG
     k_Bio_AKG = 7.0416e+03;         %
     b_Bio_AKG = 0.002;
   % OAA
     k_Bio_OAA = 46080.0;            %
     b_Bio_OAA = 0.002;

%% Kinetic parameters for Cya, degradation of cAMP, and association/dissociation ofTF
% Cya
     Cya_Vmax = 0.993000;                 % [É mol/gDCW.s]
     Cya_KEIIA = 0.0017;                  % [gProt/gDCW]
% cAMPdegr
     cAMPdegr_Vmax = 1.0;                 % [É mol/gDCW.s]
     cAMPdegr_KcAMP = 0.100000;           % [É mol/gDCW]
% cAMP-Crp
     TF_cAMPCrp_scale = 100000000.000000; % [/s]
     TF_cAMPCrp_Kcamp = 0.895000;         % [É mol/gDCW]
     TF_cAMPCrp_n = 1.000000;             %
% Cra-FBP
     TF_CraFBP_scale = 100.000000;        % [/s]
     TF_CraFBP_Kfbp = 1.360000;           % [É mol/gDCW]
     TF_CraFBP_n = 2.000000;              %
% PdhR-PYR
     TF_PdhRPYR_scale = 100.000000;       % [/s]
     TF_PdhRPYR_Kpyr = 0.164000;          % [É mol/gDCW]
     TF_PdhRPYR_n = 1.000000;             %

%% Kinetic parameters for electron transport chain
% Nuo
     v_Nuo_max = 708.3333;  % [É mol/gDCW.s]
     Km_Nuo_Q = 0.34;       % [É mol/gDCW]
     Km_Nuo_e2H2 = 0.02e+3; % [É mol/gDCW]
% Ndh
     v_Ndh_max = 11.1111;   % [É mol/gDCW.s]
     Km_Ndh_Q = 0.34;       % [É mol/gDCW]
     Km_Ndh_e2H2 = 0.02e+3; % [É mol/gDCW]
% Cyo
     v_Cyo_max = 800.0;     % [É mol/gDCW.s]
     Km_Cyo_O2 = 2.0e-4;    % [mM]
     Km_Cyo_QH2 = 25.0;     % [É mol/gDCW]
     k_O2 = 0.001;
% Cyd
     v_Cyd_max = 31.5582;   % [É mol/gDCW.s]
     Km_Cyd_O2 = 2.4e-5;    % [mM]
     Km_Cyd_QH2 = 10.0;     % [É mol/gDCW]
% O2
    O2_saturation = 0.214;  % [mM]
% Synthesis of QH2
    p_syn_QH2 = 2.16;       % [É mol/gDCW.s]

%% Effect of TF on gene expression
% cAMP-Crp
    K_CrpcAMP = 0.001200;   % [gProt/gDCW]
% Cra
    K_Cra = 0.015000;       % [gProt/gDCW]
% PdhR
    K_PdhR = 0.001;         % [gProt/gDCW]
% IclR
    n_IclR = -3;            %
    K_IclR = 0.003;         % [gProt/gDCW]
% ArcA
    K_A_ArcA_Q = 0.5;       % [É mol/gDCW]
    n_ArcA = -5.0;          %
% Fnr
    K_A_FNR_O2 = 5.0e-6;    % [mM]
    n_FNR = -5.0;           %
    
%% Total concentrations for TF and cofactor
% Crp_total
    Crp_total = 7.29e-03;   % [gProt/gDCW]
% Cra_total
    Cra_total = 7.29e-03;   % [gProt/gDCW]
% PdhR_total
    PdhR_total = 7.29e-03;  % [gProt/gDCW]
% IclR
    IclR = 7.29e-03;        % [gProt/gDCW]
% ATP
    ATP = 7.5709;           % [É mol/gDCW]
% NADP
%   NADP = 0.3457;          % [É mol/gDCW]
    NADP = 0.0037;          % [É mol/gDCW]
% NADPH
%     NADPH = 0.1099;         % [É mol/gDCW]
% NA*_total
    NADH_NAD_total = 12.90; % [É mol/gDCW]
% CoA
    CoA = 0.8865;           % [É mol/gDCW]
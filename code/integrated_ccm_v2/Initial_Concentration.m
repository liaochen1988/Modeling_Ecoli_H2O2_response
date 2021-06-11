function ic = Initial_Concentration()

% Fermentative variables
   OD_ini = 0.030000;     % [-]
   GLC_ex_ini = 10.0;     % [g/l]
   ACT_ex_ini = 0.000000; % [g/l]
   LAC_ini = 0.0;         % [g/l]
   FOR_ini = 0.0;         % [g/l]
   ETH_ini = 0.0;         % [g/l]
% Intracellular metabolites
 % Glycolysis
   GLCin_ini = 0.001000;  % [��mol/gDCW]
   G6P_ini = 1.908141;    % [��mol/gDCW]
   FBP_ini = 6.575042;    % [��mol/gDCW]
   PG3_ini = 5.720977;    % [��mol/gDCW]
   PEP_ini = 0.210456;    % [��mol/gDCW]
   PYR_ini = 0.863278;    % [��mol/gDCW]
   ACoA_ini = 0.351972;   % [��mol/gDCW]
   AcAld_ini = 0.1950;    % [��mol/gDCW]
 % PP pathway
   PG6_ini = 1.4326;      % [��mol/gDCW]
   Ribu5_ini = 0.1968;    % [��mol/gDCW]
   Rib5_ini = 0.7057;     % [��mol/gDCW]
   Xyl5_ini = 0.2447;     % [��mol/gDCW]
   S7P_ini = 0.4894;      % [��mol/gDCW]
   E4P_ini = 0.1738;      % [��mol/gDCW]
 % TCA cycle
   ICT_ini = 0.001408116; % [��mol/gDCW]
   AKG_ini = 0.191191;    % [��mol/gDCW]
   SUC_ini = 3.278779;    % [��mol/gDCW]
   MAL_ini = 3.278779;    % [��mol/gDCW]
   OAA_ini = 0.050535;    % [��mol/gDCW]
   GLX_ini = 5.70593E-09; % [��mol/gDCW]
% PTS protein
   EIIA_ini = 0.096477;   % [��mol/gDCW]
   EIIA_P_ini = 0.003523; % [��mol/gDCW]
% cAMP and transcription factor
   cAMP_ini = 0.202804;   % [��mol/gDCW]
   CrpcAMP_ini = 0.001347;% [��mol/gDCW]
   CraFBP_ini = 0.006991; % [��mol/gDCW]
   PdhRPYR_ini = 0.006126;% [��mol/gDCW]
% Electron transport chain
   QH2_ini = 1.0;         % [��mol/gDCW]
   Q_ini = 1.0;           % [��mol/gDCW]
   NADH_ini = 0.1475;     % umol/gDCW, Bennett2009NatChemBiol, 4 g/L glucose MM, 83.2 uM
   NAD_ini = 4.5213;      % umol/gDCW, Bennett2009NatChemBiol, 4 g/L glucose MM, 2550 uM
% Extracellular SUC conc.
   SUC_ex_ini = 0.0;      % [g/l]
% ATP conc.
   ATP_ini = 7.5709;      % umol/gDCW
   NADP_ini = 0.0037;     % umol/gDCW, Bennett2009NatChemBiol, 4 g/L glucose MM, 2.08 uM
   NADPH_ini = 0.2145;    % umol/gDCW, Bennett2009NatChemBiol, 4 g/L glucose MM, 121 uM
   MET_ini = 1.0;
   MetSO_ini = 0.0;
   ppGpp_ini = 1.0;
% Cellular protein (M)
   PSH_ini = 3.0;
   PSS_ini = 0.0;
% MetE
   MetE_SH_ini = 1.0;
   MetE_SSG_ini = 0.0;
   RIB_ini = 1.0;
   H2O2_ex_ini = 0.0;
   H2O2_ini = 0.0;
% OxyR conc. (M)
   OxyR_ox_ini = 0.0;
   OxyR_red_ini = 1.96e-6; % Estimated from Li 2014, MOPS minimal medium
% Catalase concentration
   KatE_ini = 5e-6;
   KatG_ini = 0.0;
   AhpCF_ini = 0.0;
   Grx_ini = 1e-5;
   Gor_ini = 0.0;
% GSH, GSSG (M)
   GSH_ini  = 0.0166;    % Bennett2009NatChemBiol, 4 g/L glucose MM
   GSSG_ini = 0.00237;   % Bennett2009NatChemBiol, 4 g/L glucose MM
% Thioredoxin (M)
   Trx_ox_ini = 0.0;
   Trx_red_ini = 1.0e-5; % AdolfsenBrynildsen2015PCB
   
%% vectorize all Init conc.

ic = [OD_ini GLC_ex_ini ACT_ex_ini G6P_ini FBP_ini PG3_ini PEP_ini PYR_ini ACoA_ini ICT_ini AKG_ini MAL_ini OAA_ini GLX_ini PG6_ini Ribu5_ini Xyl5_ini Rib5_ini S7P_ini E4P_ini GLCin_ini EIIA_ini EIIA_P_ini cAMP_ini CrpcAMP_ini CraFBP_ini PdhRPYR_ini SUC_ini LAC_ini AcAld_ini ETH_ini FOR_ini QH2_ini Q_ini NADH_ini SUC_ex_ini ...
      ATP_ini NADPH_ini MET_ini MetSO_ini ppGpp_ini PSH_ini PSS_ini MetE_SH_ini MetE_SSG_ini RIB_ini H2O2_ex_ini H2O2_ini ...
      OxyR_ox_ini OxyR_red_ini KatE_ini KatG_ini AhpCF_ini Grx_ini Gor_ini GSH_ini GSSG_ini Trx_ox_ini Trx_red_ini NADP_ini NAD_ini];


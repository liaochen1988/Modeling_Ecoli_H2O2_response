function ic = Initial_Concentration(GLC_ex_ini)

% Fermentative variables
   OD_ini = 0.030000;     % [-]
   %GLC_ex_ini = 10.0;     % [g/l]
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
   NADH_ini = 6.64;       % [��mol/gDCW]
% Extracellular SUC conc.
   SUC_ex_ini = 0.0;      % [g/l]
   
   
%% vectorize all Init conc.

ic = [OD_ini GLC_ex_ini ACT_ex_ini G6P_ini FBP_ini PG3_ini PEP_ini PYR_ini ACoA_ini ICT_ini AKG_ini MAL_ini OAA_ini GLX_ini PG6_ini Ribu5_ini Xyl5_ini Rib5_ini S7P_ini E4P_ini GLCin_ini EIIA_ini EIIA_P_ini cAMP_ini CrpcAMP_ini CraFBP_ini PdhRPYR_ini SUC_ini LAC_ini AcAld_ini ETH_ini FOR_ini QH2_ini Q_ini NADH_ini SUC_ex_ini];


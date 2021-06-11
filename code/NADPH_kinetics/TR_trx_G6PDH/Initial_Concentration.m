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
   GLCin_ini = 0.001000;  % [É mol/gDCW]
   G6P_ini = 1.908141;    % [É mol/gDCW]
   FBP_ini = 6.575042;    % [É mol/gDCW]
   PG3_ini = 5.720977;    % [É mol/gDCW]
   PEP_ini = 0.210456;    % [É mol/gDCW]
   PYR_ini = 0.863278;    % [É mol/gDCW]
   ACoA_ini = 0.351972;   % [É mol/gDCW]
   AcAld_ini = 0.1950;    % [É mol/gDCW]
 % PP pathway
   PG6_ini = 1.4326;      % [É mol/gDCW]
   Ribu5_ini = 0.1968;    % [É mol/gDCW]
   Rib5_ini = 0.7057;     % [É mol/gDCW]
   Xyl5_ini = 0.2447;     % [É mol/gDCW]
   S7P_ini = 0.4894;      % [É mol/gDCW]
   E4P_ini = 0.1738;      % [É mol/gDCW]
 % TCA cycle
   ICT_ini = 0.001408116; % [É mol/gDCW]
   AKG_ini = 0.191191;    % [É mol/gDCW]
   SUC_ini = 3.278779;    % [É mol/gDCW]
   MAL_ini = 3.278779;    % [É mol/gDCW]
   OAA_ini = 0.050535;    % [É mol/gDCW]
   GLX_ini = 5.70593E-09; % [É mol/gDCW]
% PTS protein
   EIIA_ini = 0.096477;   % [É mol/gDCW]
   EIIA_P_ini = 0.003523; % [É mol/gDCW]
% cAMP and transcription factor
   cAMP_ini = 0.202804;   % [É mol/gDCW]
   CrpcAMP_ini = 0.001347;% [É mol/gDCW]
   CraFBP_ini = 0.006991; % [É mol/gDCW]
   PdhRPYR_ini = 0.006126;% [É mol/gDCW]
% Electron transport chain
   QH2_ini = 1.0;         % [É mol/gDCW]
   Q_ini = 1.0;           % [É mol/gDCW]
   NADH_ini = 6.64;       % [É mol/gDCW]
% Extracellular SUC conc.
   SUC_ex_ini = 0.0;      % [g/l]
% NADPH
   NADPH_ini = 123/564;   % [Emol/gDCW]
% Trx
   Trx_ox_ini = 1.4e-6;   % M
   Trx_red_ini = 8.6e-6;  % M
   
%% vectorize all Init conc.

ic = [OD_ini GLC_ex_ini ACT_ex_ini G6P_ini FBP_ini PG3_ini PEP_ini PYR_ini ACoA_ini ICT_ini AKG_ini MAL_ini OAA_ini GLX_ini PG6_ini Ribu5_ini Xyl5_ini Rib5_ini S7P_ini E4P_ini GLCin_ini EIIA_ini EIIA_P_ini cAMP_ini CrpcAMP_ini CraFBP_ini PdhRPYR_ini SUC_ini LAC_ini AcAld_ini ETH_ini FOR_ini QH2_ini Q_ini NADH_ini SUC_ex_ini NADPH_ini Trx_ox_ini Trx_red_ini];


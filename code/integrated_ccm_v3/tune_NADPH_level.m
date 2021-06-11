clear all;
clc;

%% Simulatioon options

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose (g/L)
glucose = 10 * 1e-3 * 180.156;

%   scaled dissolved oxygen level
a = 40; % aerobic condition

Num_of_State_Variable = 61;

V_env = 0.02; % L

%% Simulation
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);
IC = Initial_Concentration();
IC(2) = glucose;

% wild type
katG_mutant = false;
katE_mutant = false;
const_ext_vars = ones(1,8);
[t,x] = ode15s(@Kinetic_model,[-10000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true,const_ext_vars);
res = Kinetic_model(t(end),x(end,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false,const_ext_vars);
res(1) * 3600 % desired growth rate (0.9 per hour)

%% Print
uc_gDCW_L = 564; % [gDCW/l]
NADH = x(end, 35) * uc_gDCW_L;
ATP = x(end, 37) * uc_gDCW_L * 1e-3;
NADPH = x(end, 38) * uc_gDCW_L;
NADP = x(end, 60) * uc_gDCW_L;
GSH = x(end, 56) * 1e3;
GSSG = x(end, 57) * 1e3;
fprintf('ATP: %f mM\n', ATP); % ATP 9.63, [8.13, 11.4]
fprintf('NADH: %f uM\n', NADH); % NADH [54.5, 127]
fprintf('NADP: %f uM\n', NADP); % NADP 2.08, [0.14, 31.1]
fprintf('NADPH: %f uM\n', NADPH); % NADPH 121, [110, 134]
fprintf('GSH: %f mM\n', GSH); % GSH 16.6, [15.3, 17.9]
fprintf('GSSG: %f mM\n', GSSG); % GSSG 2.37, [1.94, 2.90]
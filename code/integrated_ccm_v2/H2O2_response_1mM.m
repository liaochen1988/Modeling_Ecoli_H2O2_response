clear all;
clc;

%% Experimental data
t_1mM = [0, 5, 15, 30, 60]; % sec
G6P = [2.085	1.689	1.849	1.670	1.561]; % mM
F6P = [0.294	0.254	0.313	0.256	0.265]; % mM
FBP = [1.132	1.237	1.522	1.579	1.134]; % mM
DHAP = [0.417	0.426	0.470	0.468	0.446]; % mM
PEP = [0.042	0.028	0.017	0.023	0.016]; % mM
AcCoA = [1.000	1.272	1.453	1.526	1.642]; % Relative
Aconitate = [1.000	1.385	2.233	3.662	6.086]; % Relative
AKG = [0.181	0.156	0.108	0.088	0.132]; % mM
Malate = [0.853	0.573	0.389	0.269	0.250]; % mM
PG6 = [0.137	0.319	0.262	0.196	0.168]; % mM
Ru5P = [0.444	0.513	0.533	0.431	0.442]; % mM
XU5P = [0.335	0.466	0.411	0.403	0.293]; % mM
R5P = [0.462	0.608	0.565	0.490	0.465]; % mM
S7P = [0.096	0.138	0.137	0.122	0.120]; % mM
Succinate = [0.149	0.141	0.116	0.093	0.087]; % mM
NADH = [1.000	0.668	0.846	0.790	0.898]; % Relative
NADPH = [1.000	1.163	0.645	0.697	0.568]; % Relative

% P5P = [0.980, 1.296, 1.221, 1.049, 0.999];
% Succoa = [1.000	1.253	0.621	0.483	0.539]; % Relative
% ADP = [0.382	0.339	0.413	0.309	0.274]; % mM
% AMP = [0.395	0.464	0.388	0.303	0.323]; % mM
% NAD = [1.000	2.471	2.492	1.914	2.391]; % Relative
% NADP = [1.000	1.063	0.901	0.916	0.907]; % Relative
% XPG = [1.698	1.800	1.494	1.350	1.264]; % mM

%% Simulatioon options

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose (g/L)
glucose = 2;

%   scaled dissolved oxygen level
a = 40; % aerobic condition

Num_of_State_Variable = 59;

tspan = [0:0.1:60];
H2O2_input = 1e-3; % M
OD = 0.5;
N_cell_per_OD = 1e9; % cells/mL/OD
V_env = 0.02; % L

%% Initial condition
% tolerance cannot be too tight
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);
IC = Initial_Concentration();
IC(2) = glucose;

% wild type
katG_mutant = false;
katE_mutant = false;
[~,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
IC_wt = x(end,:);

%% simulation
options=odeset('RelTol',1e-10,'AbsTol',1e-10, 'NonNegative',[1:Num_of_State_Variable]);
IC_wt(1) = OD * N_cell_per_OD * V_env * 1e3;
IC_wt(47) = H2O2_input;
[t_wt,x_wt] = ode15s(@Kinetic_model,tspan,IC_wt,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
mu = zeros(length(t_wt),1);
for i=1:length(t_wt)
    res = Kinetic_model(t_wt(i),x_wt(i,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false);
    mu(i) = res(1);
end

%% plot 
figure();

subplot(2,3,1);
hold on;
plot(t_wt, x_wt(:,48) * 1e6,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('[H_2O_2]_{in} (\muM)');

subplot(2,3,2);
hold on;
plot(t_wt, x_wt(:,47) * 1e6,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('[H_2O_2]_{ex} (\muM)');

subplot(2,3,3);
hold on;
plot(t_wt, mu * 3600,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('Growth rate (h^{-1})');

subplot(2,3,4);
hold on;
plot(t_wt, x_wt(:,56) * 1e3,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('GSH (mM)');

subplot(2,3,5);
hold on;
plot(t_wt, x_wt(:,57) * 1e3, 'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('GSSG (mM)');

subplot(2,3,6);
hold on;
plot(t_wt, x_wt(:,55) * 1e6,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('Gor (\muM)');









figure();

uc_gDCW_L = 564;         % [gDCW/l]

% G6P
subplot(4,4,1);
hold on;
plot(tspan, x_wt(:,4) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, G6P, 'ro');
plot(t_1mM, F6P, 'bo');
axis square;
box on;
xlabel('Time (s)');
ylabel('G6P (mM)');

% FBP
subplot(4,4,2);
hold on;
plot(tspan, x_wt(:,5) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, FBP, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('FBP (mM)');

% GAP/DHAP
subplot(4,4,3);
hold on;
plot(tspan, x_wt(:,6) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, DHAP, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('PEP (mM)');

% PEP
subplot(4,4,4);
hold on;
plot(tspan, x_wt(:,7) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, PEP, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('PEP (mM)');

% ACoA
subplot(4,4,5);
hold on;
plot(tspan, x_wt(:,9) ./ x_wt(1,9),'k-');
plot(t_1mM, AcCoA, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('AcCoA');

% Aconitate
subplot(4,4,6);
hold on;
plot(tspan, x_wt(:,10) ./ x_wt(1,10),'k-');
plot(t_1mM, Aconitate, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Aconitate');

% Aconitate
subplot(4,4,7);
hold on;
plot(tspan, x_wt(:,11) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, AKG, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('AKG (mM)');

% Malate
subplot(4,4,8);
hold on;
plot(tspan, x_wt(:,12) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, Malate, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Malate (mM)');

% 6PG
subplot(4,4,9);
hold on;
plot(tspan, x_wt(:,15) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, PG6, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('6PG (mM)');

% Ru5P
subplot(4,4,10);
hold on;
plot(tspan, x_wt(:,16) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, Ru5P, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Ru5P (mM)');

% Xu5P
subplot(4,4,11);
hold on;
plot(tspan, x_wt(:,17) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, XU5P, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('XU5P (mM)');

% R5P
subplot(4,4,12);
hold on;
plot(tspan, x_wt(:,18) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, R5P, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('R5P (mM)');

% S7P
subplot(4,4,13);
hold on;
plot(tspan, x_wt(:,19) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, S7P, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('S7P (mM)');

% Succinate
subplot(4,4,14);
hold on;
plot(tspan, x_wt(:,28) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, Succinate, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Succinate (mM)');

% NADH
subplot(4,4,15);
hold on;
plot(tspan, x_wt(:,35) ./ x_wt(1,35),'k-');
plot(t_1mM, NADH, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('NADH');

% NADPH
subplot(4,4,16);
hold on;
plot(tspan, x_wt(:,38) ./ x_wt(1,38),'k-');
plot(t_1mM, NADPH, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('NADPH');
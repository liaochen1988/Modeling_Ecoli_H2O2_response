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

Num_of_State_Variable = 62;

tspan = [0:0.1:60];
H2O2_input = 1e-3; % M
OD = 0.5;
N_cell_per_OD = 1e9; % cells/mL/OD
V_env = 0.02; % L

%% Initial condition
% tolerance cannot be too tight
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);
IC = Initial_Concentration();
IC(1) = OD * N_cell_per_OD * V_env * 1e3;
IC(2) = glucose;

% wild type
katG_mutant = false;
katE_mutant = false;
const_ext_vars = zeros(1,10);
const_ext_vars([1,2]) = 1; % N and GLC_ex
[t,x] = ode15s(@Kinetic_model,[-10000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true,const_ext_vars);
res = Kinetic_model(t(end),x(end,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false,const_ext_vars);
res(1) * 3600
res(3)
x(end,54) % Grx
x(end,56) % GSH
x(end,57) % GSSG
x(end,38) % nadph
x(end,44)/(x(end,45)+x(end,44))
IC_wt = x(end,:);

%% simulation
options=odeset('RelTol',1e-10,'AbsTol',1e-10, 'NonNegative',[1:Num_of_State_Variable]);
IC_wt(47) = H2O2_input;
IC_wt(38) = 121/564; % 0.2145
const_ext_vars = zeros(1,10);
[t_wt,x_wt] = ode15s(@Kinetic_model,tspan,IC_wt,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true,const_ext_vars);
mu = zeros(length(t_wt),1);
v_ICDH = zeros(length(t_wt),1);
v_CS = zeros(length(t_wt),1);
v_Icl = zeros(length(t_wt),1);
v_GSH_red2 = zeros(length(t_wt),1);
for i=1:length(t_wt)
    res = Kinetic_model(t_wt(i),x_wt(i,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false,const_ext_vars);
    mu(i) = res(1);
    v_ICDH(i) = res(4);
    v_CS(i) = res(5);
    v_Icl(i) = res(6);
    v_GSH_red2(i) = res(7);
end

%% plot 
figure();

subplot(1,2,1)
hold on;
plot(t_wt/60, -(-252 - 8.314 * 310 / 2 / 96485 * log(x_wt(:,56).^2 ./ x_wt(:,57)) * 1000), 'k-');
axis square;
box on;
xlabel('Time (min)');
ylabel('Redox potential');
xlim([-1,6]);

subplot(1,2,2);
hold on;
plot(t_wt/60, (3-x_wt(:,42))/3, 'k-');
axis square;
box on;
xlabel('Time (min)');
ylabel('Oxidized protein (%)');
xlim([-1,6]);

figure();

subplot(3,4,1);
hold on;
plot(t_wt, x_wt(:,48) * 1e6,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('[H_2O_2]_{in} (\muM)');

subplot(3,4,2);
hold on;
plot(t_wt, x_wt(:,47) * 1e6,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('[H_2O_2]_{ex} (\muM)');

subplot(3,4,3);
hold on;
plot(t_wt, mu * 3600,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('Growth rate (h^{-1})');

subplot(3,4,4);
hold on;
plot(t_wt, x_wt(:,56) * 1e3,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('GSH (mM)');

subplot(3,4,5);
hold on;
plot(t_wt, x_wt(:,57) * 1e3, 'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('GSSG (mM)');

subplot(3,4,6);
hold on;
plot(t_wt, x_wt(:,55) * 1e6,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('Gor (\muM)');

subplot(3,4,7);
hold on;
plot(t_wt, x_wt(:,62),'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('ICDH activity');

subplot(3,4,8);
hold on;
plot(t_wt, v_ICDH,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('ICDH reaction flux');

subplot(3,4,9);
hold on;
plot(t_wt, v_CS,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('CS reaction flux');

subplot(3,4,10);
hold on;
plot(t_wt, v_Icl,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('Icl reaction flux');

subplot(3,4,11);
hold on;
plot(t_wt, mu .* x_wt(:,10),'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('mu * ICT');

subplot(3,4,12);
hold on;
plot(t_wt, v_GSH_red2,'k-');
axis square;
box on;
xlabel('Time (s)');
ylabel('v_GSH_red2');

figure();

uc_gDCW_L = 564;         % [gDCW/l]

% G6P
subplot(4,4,1);
hold on;
% plot(tspan, x_wt(:,4) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, G6P, 'ro');
plot(tspan, x_wt(:,4) / x_wt(1,4),'k-');
plot(t_1mM, G6P / G6P(1), 'ro');
%plot(t_1mM, F6P, 'bo');
axis square;
box on;
xlabel('Time (s)');
ylabel('G6P (mM)');

% FBP
subplot(4,4,2);
hold on;
% plot(tspan, x_wt(:,5) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, FBP, 'ko');
plot(tspan, x_wt(:,5) / x_wt(1,5),'k-');
plot(t_1mM, FBP / FBP(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('FBP (mM)');

% GAP/DHAP
subplot(4,4,3);
hold on;
% plot(tspan, x_wt(:,6) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, DHAP, 'ko');
plot(tspan, x_wt(:,6) / x_wt(1,6),'k-');
plot(t_1mM, DHAP / DHAP(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('PEP (mM)');

% PEP
subplot(4,4,4);
hold on;
% plot(tspan, x_wt(:,7) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, PEP, 'ko');
plot(tspan, x_wt(:,7) / x_wt(1,7),'k-');
plot(t_1mM, PEP / PEP(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('PEP (mM)');

% ACoA
subplot(4,4,5);
hold on;
plot(tspan, x_wt(:,9) / x_wt(1,9),'k-');
plot(t_1mM, AcCoA, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('AcCoA');

% Aconitate
subplot(4,4,6);
hold on;
plot(tspan, x_wt(:,10)/ x_wt(1,10),'k-');
plot(t_1mM, Aconitate, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Aconitate');

% Aconitate
subplot(4,4,7);
hold on;
% plot(tspan, x_wt(:,11) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, AKG, 'ko');
plot(tspan, x_wt(:,11) / x_wt(1,11),'k-');
plot(t_1mM, AKG / AKG(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('AKG (mM)');

% Malate
subplot(4,4,8);
hold on;
% plot(tspan, x_wt(:,12) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, Malate, 'ko');
plot(tspan, x_wt(:,12) / x_wt(1,12),'k-');
plot(t_1mM, Malate / Malate(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Malate (mM)');

% 6PG
subplot(4,4,9);
hold on;
% plot(tspan, x_wt(:,15) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, PG6, 'ko');
plot(tspan, x_wt(:,15) / x_wt(1,15),'k-');
plot(t_1mM, PG6 / PG6(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('6PG (mM)');

% Ru5P
subplot(4,4,10);
hold on;
% plot(tspan, x_wt(:,16) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, Ru5P, 'ko');
plot(tspan, x_wt(:,16) / x_wt(1,16), 'k-');
plot(t_1mM, Ru5P / Ru5P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Ru5P (mM)');

% Xu5P
subplot(4,4,11);
hold on;
% plot(tspan, x_wt(:,17) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, XU5P, 'ko');
plot(tspan, x_wt(:,17) / x_wt(1,17), 'k-');
plot(t_1mM, XU5P / XU5P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('XU5P (mM)');

% R5P
subplot(4,4,12);
hold on;
% plot(tspan, x_wt(:,18) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, R5P, 'ko');
plot(tspan, x_wt(:,18) / x_wt(1,18),'k-');
plot(t_1mM, R5P / R5P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('R5P (mM)');

% S7P
subplot(4,4,13);
hold on;
% plot(tspan, x_wt(:,19) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, S7P, 'ko');
plot(tspan, x_wt(:,19) / x_wt(1,19),'k-');
plot(t_1mM, S7P / S7P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('S7P (mM)');

% Succinate
subplot(4,4,14);
hold on;
% plot(tspan, x_wt(:,28) * uc_gDCW_L * 1e-3,'k-');
% plot(t_1mM, Succinate, 'ko');
plot(tspan, x_wt(:,28) / x_wt(1,28), 'k-');
plot(t_1mM, Succinate / Succinate(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Succinate (mM)');

% NADH
subplot(4,4,15);
hold on;
plot(tspan, x_wt(:,35) / x_wt(1,35), 'k-');
plot(t_1mM, NADH, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('NADH');

% NADPH
subplot(4,4,16);
hold on;
plot(tspan, x_wt(:,38) / x_wt(1,38), 'k-');
plot(t_1mM, NADPH, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('NADPH');
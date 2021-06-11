clear all;
clc;

% whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

% initial external glucose (g/L)
glucose = 2;

% scaled dissolved oxygen level
a = 40; % aerobic condition

Num_of_State_Variable = 39;
t = [0:0.1:60];

options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

% Run NADPH response
TrxR = 10 * 1e-6; % M
kd_nadph = 1;
k_Trxox_H2O2 = 2e4; % M/s

IC = Initial_Concentration();
IC(2) = glucose;

[~,x] = ode15s(@Kinetic_model,[-10000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,true,kd_nadph,TrxR,0,k_Trxox_H2O2);
IC = x(end,:);
IC(37) * 564
[t,x] = ode15s(@Kinetic_model,[0.0 60.0],IC,options,arcA_mutant,fnr_mutant,a,true,kd_nadph,TrxR,1e-3,k_Trxox_H2O2);
v_G6PDH = zeros(1,length(t));
v_PGDH = zeros(1,length(t));
mu = zeros(1,length(t));
for i=1:length(t)
    res = Kinetic_model(t(i),x(i,:),arcA_mutant,fnr_mutant,a,false,kd_nadph,TrxR,1e-3,k_Trxox_H2O2);
    v_G6PDH(i) = res(1);
    v_PGDH(i) = res(2);
    mu(i) = res(3);
end

%% Experimental data
t_1mM = [0, 5, 15, 30, 60]; % sec
G6P = [2.085	1.689	1.849	1.670	1.561]; % mM
G6P_sd= [0.963	0.589	0.752	0.733	0.682];
F6P = [0.294	0.254	0.313	0.256	0.265]; % mM
FBP = [1.132	1.237	1.522	1.579	1.134]; % mM
DHAP = [0.417	0.426	0.470	0.468	0.446]; % mM
PEP = [0.042	0.028	0.017	0.023	0.016]; % mM
AcCoA = [1.000	1.272	1.453	1.526	1.642]; % Relative
Aconitate = [1.000	1.385	2.233	3.662	6.086]; % Relative
AKG = [0.181	0.156	0.108	0.088	0.132]; % mM
Malate = [0.853	0.573	0.389	0.269	0.250]; % mM
PG6 = [0.137	0.319	0.262	0.196	0.168]; % mM
PG6_sd = [0.030	0.105	0.078	0.074	0.087];
Ru5P = [0.444	0.513	0.533	0.431	0.442]; % mM
XU5P = [0.335	0.466	0.411	0.403	0.293]; % mM
R5P = [0.462	0.608	0.565	0.490	0.465]; % mM
S7P = [0.096	0.138	0.137	0.122	0.120]; % mM
Succinate = [0.149	0.141	0.116	0.093	0.087]; % mM
NADH = [1.000	0.668	0.846	0.790	0.898]; % Relative
NADPH = [1.000	1.163	0.645	0.697	0.568]; % Relative
NADPH_sd = [0.108	0.601	0.246	0.248	0.192];

figure();

uc_gDCW_L = 564;         % [gDCW/l]

% G6P
subplot(4,4,1);
hold on;
plot(t, x(:,4) * uc_gDCW_L * 1e-3,'k-');
errorbar(t_1mM, G6P, G6P_sd,'ko');
% plot(t, x(:,4) / x(1,4),'k-');
% plot(t_1mM, G6P / G6P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('G6P (mM)');
axis([-6,60,0,4]);

% 6PG
subplot(4,4,2);
hold on;
plot(t, x(:,15) * uc_gDCW_L * 1e-3,'k-');
errorbar(t_1mM, PG6, PG6_sd, 'ko');
% plot(t, x(:,15) / x(1,15),'k-');
% plot(t_1mM, PG6 / PG6(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('6PG (mM)');
%axis([-6,60,0,0.5]);

% FBP
subplot(4,4,3);
hold on;
plot(t, x(:,5) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, FBP, 'ko');
% plot(t, x(:,5) / x(1,5),'k-');
% plot(t_1mM, FBP / FBP(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('FBP (mM)');

% GAP/DHAP
subplot(4,4,4);
hold on;
plot(t, x(:,6) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, DHAP, 'ko');
% plot(t, x(:,6) / x(1,6),'k-');
% plot(t_1mM, DHAP / DHAP(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('PEP (mM)');

% PEP
subplot(4,4,5);
hold on;
plot(t, x(:,7) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, PEP, 'ko');
% plot(t, x(:,7) / x(1,7),'k-');
% plot(t_1mM, PEP / PEP(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('PEP (mM)');

% ACoA
subplot(4,4,6);
hold on;
plot(t, x(:,9) / x(1,9),'k-');
plot(t_1mM, AcCoA, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('AcCoA');

% Aconitate
subplot(4,4,7);
hold on;
plot(t, x(:,10)/ x(1,10),'k-');
plot(t_1mM, Aconitate, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Aconitate');

% Aconitate
subplot(4,4,8);
hold on;
plot(t, x(:,11) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, AKG, 'ko');
% plot(t, x(:,11) / x(1,11),'k-');
% plot(t_1mM, AKG / AKG(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('AKG (mM)');

% Malate
subplot(4,4,9);
hold on;
plot(t, x(:,12) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, Malate, 'ko');
% plot(t, x(:,12) / x(1,12),'k-');
% plot(t_1mM, Malate / Malate(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Malate (mM)')

% Ru5P
subplot(4,4,10);
hold on;
plot(t, x(:,16) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, Ru5P, 'ko');
% plot(t, x(:,16) / x(1,16), 'k-');
% plot(t_1mM, Ru5P / Ru5P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Ru5P (mM)');

% Xu5P
subplot(4,4,11);
hold on;
plot(t, x(:,17) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, XU5P, 'ko');
% plot(t, x(:,17) / x(1,17), 'k-');
% plot(t_1mM, XU5P / XU5P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('XU5P (mM)');

% R5P
subplot(4,4,12);
hold on;
% % plot(t, x(:,18) * uc_gDCW_L * 1e-3,'k-');
% % plot(t_1mM, R5P, 'ko');
plot(t, x(:,18) / x(1,18),'k-');
plot(t_1mM, R5P / R5P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('R5P (mM)');

% S7P
subplot(4,4,13);
hold on;
plot(t, x(:,19) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, S7P, 'ko');
% plot(t, x(:,19) / x(1,19),'k-');
% plot(t_1mM, S7P / S7P(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('S7P (mM)');

% Succinate
subplot(4,4,14);
hold on;
plot(t, x(:,28) * uc_gDCW_L * 1e-3,'k-');
plot(t_1mM, Succinate, 'ko');
% plot(t, x(:,28) / x(1,28), 'k-');
% plot(t_1mM, Succinate / Succinate(1), 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('Succinate (mM)');

% NADH
subplot(4,4,15);
hold on;
plot(t, x(:,35) / x(1,35), 'k-');
plot(t_1mM, NADH, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('NADH');

% NADPH
subplot(4,4,16);
hold on;
plot(t, x(:,37) / x(1,37), 'k-');
errorbar(t_1mM, NADPH, NADPH_sd, 'ko');
axis square;
box on;
xlabel('Time (s)');
ylabel('NADPH');
axis([-6,60,0,2]);
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
tspan = [0:0.1:60];

options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

% Run NADPH response
TrxR = 8.05 * 1e-6; % M
kd_nadph = 3.05;
k_Trxox_H2O2 = 10.^[0:0.1:8];
resnorm = zeros(length(k_Trxox_H2O2),1);

% data
t_1mM = [0, 5, 15, 30, 60]; % sec
NADPH_1mM = [1.000	1.163	0.645	0.697	0.568]; % Relative
NADPH_1mM_sd = [0.108	0.601	0.246	0.248	0.192];

for i=1:length(k_Trxox_H2O2)
    i
    IC = Initial_Concentration();
    IC(2) = glucose;
    
    [~,x] = ode15s(@Kinetic_model,[-10000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,true,kd_nadph,TrxR,0,k_Trxox_H2O2(i));
    IC = x(end,:);
    [t,x] = ode15s(@Kinetic_model,[0.0 60.0],IC,options,arcA_mutant,fnr_mutant,a,true,kd_nadph,TrxR,1e-3,k_Trxox_H2O2(i));
    
    resnorm(i) = sum(((pchip(t, x(:,37)/x(1,37), t_1mM) - NADPH_1mM)./NADPH_1mM_sd).^2);
end

% best fit: k_H2O2 = 5e4

figure();
hold on;

plot(k_Trxox_H2O2, resnorm);
axis square;
box on;
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
k_Trxox_H2O2 = [0, 1e5, 1e6, 1e7];
cls = hsv(length(k_Trxox_H2O2));

figure();

for i=1:length(k_Trxox_H2O2)
    IC = Initial_Concentration();
    IC(2) = glucose;
    
    [t,x] = ode15s(@Kinetic_model,[-10000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,true,kd_nadph,TrxR,0,k_Trxox_H2O2(i));
    IC = x(end,:);
    IC(37) * 564
    [t,x] = ode15s(@Kinetic_model,[0.0 60.0],IC,options,arcA_mutant,fnr_mutant,a,true,kd_nadph,TrxR,1e-3,k_Trxox_H2O2(i));
    
    subplot(1,3,1);
    hold on;
    plot(t,x(:,38) * 1e6,'k-','Color',cls(i,:));
    
    subplot(1,3,2);
    hold on;
    plot(t,x(:,37) * 564,'k-','Color',cls(i,:));
    
    subplot(1,3,3);
    hold on;
    plot(t,x(:,37)/x(1,37),'k-','Color',cls(i,:));
end

subplot(1,3,1);
hold on;
axis([-1,60,0,10]);
xlabel('Time (s)');
ylabel('Trx ox (uM)');
axis square;
box on;

subplot(1,3,2);
hold on;
axis([-1,60,0,150]);
plot([-1,60],[123,123],'k--');
xlabel('Time (s)');
ylabel('NADPH (uM)');
axis square;
box on;

subplot(1,3,3);
hold on;

axis([-1,60,0,1.8]);
xlabel('Time (s)');
ylabel('Relative NADPH');

t_1mM = [0, 5, 15, 30, 60]; % sec
NADPH_1mM = [1.000	1.163	0.645	0.697	0.568]; % Relative
NADPH_1mM_sd = [0.108	0.601	0.246	0.248	0.192];
errorbar(t_1mM, NADPH_1mM, NADPH_1mM_sd, 'ko');
axis square;
box on;
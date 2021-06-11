clear all;
clc;

%% Simulate in vitro dynamics of OxyR
% Fig. 3 in Aslund1999PNAS

kon_oxyr_h2o2 = 9.9863; % 1/s
n_oxyr = 1.3614;
Km_oxyr = 41.3179; % uM
kon_oxyr_grx = 1.4910e7; % 1/s/M^5
koff_oxyr_grx = 1.5218e15; % 1/s/M^9

oxyr_total = 1; % uM
GSH = 0.025; % M
GSSG = 1e-4; % M
Grx = 10*1e-6; % M

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,3),'NonNegative',1:3);

% before shift
x0 = [0,0,oxyr_total];
tspan = [0,1e10];
[t_bs,x_bs] = ode15s(@oxyr_kinetic_model_in_vitro,tspan,x0,options,...
    kon_oxyr_h2o2,n_oxyr,Km_oxyr,kon_oxyr_grx,koff_oxyr_grx,GSH,GSSG,Grx);

% after shift
x0 = x_bs(end,:);
x0(1) = 2; % uM
tspan = [0,10000];
[t_as,x_as] = ode15s(@oxyr_kinetic_model_in_vitro,tspan,x0,options,...
    kon_oxyr_h2o2,n_oxyr,Km_oxyr,kon_oxyr_grx,koff_oxyr_grx,GSH,GSSG,Grx);

%figure();

subplot(2,1,1);
hold on;
plot([t_bs-max(t_bs);t_as]/60,[x_bs(:,1);x_as(:,1)]);
box on;
xlabel('Time (min)');
ylabel('[H_2O_2] (\muM)');
axis([-1,60,0,2]);

subplot(2,1,2);
hold on;
plot([t_bs-max(t_bs);t_as]/60,[x_bs(:,2);x_as(:,2)]/oxyr_total);
box on;
xlabel('Time (min)');
ylabel('Oxidized OxyR fraction');
axis([-1,60,-0.1,1.1]);
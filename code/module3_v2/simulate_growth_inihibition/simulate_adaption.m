clear all;
clc;

%% parameters

P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20; % mol/s, corresponding to 14 uM/s
Km_ahp = 1.4138e-07; % M
Vc = 2*1e-15; % L
Ve = 0.2e-3;  % L
k0_cat = 4.188e6; % 1/s/M
k0_ahp = 52.4; % 1/s
n_oxyr = 1.3614;
Km_oxyr = 41.3179e-6; % M
k_oxyr_ox = 9.9863; % 1/s
k_oxyr_red = 0.0023; % 1/s
alpha_cat = 3.62e-24; % mol/s
alpha_ahp = 1.19e-22; % mol/s
lambda = 1.93/3600; % 1/s
N = 0;%0.43*5e8; 

%% simulation
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,4),'NonNegative',1:4);

% before shift
x0 = [0,0,0,0];
tspan = [0,1e10];
[t_bs,x_bs] = ode15s(@model_func,tspan,x0,options,N,...
    k_met,PA,k0_cat,k0_ahp,Km_ahp,Vc,Ve,alpha_cat,alpha_ahp,lambda,n_oxyr,k_oxyr_red,k_oxyr_ox,Km_oxyr);
f_bs = x_bs(:,2).^n_oxyr./(x_bs(:,2).^n_oxyr*(1+k_oxyr_red/k_oxyr_ox)+k_oxyr_red/k_oxyr_ox*Km_oxyr^n_oxyr);

% after shift
x0 = x_bs(end,:);
x0(1) = 1e-5; % M
tspan = [0,1e6];
[t_as,x_as] = ode15s(@model_func,tspan,x0,options,N,...
    k_met,PA,k0_cat,k0_ahp,Km_ahp,Vc,Ve,alpha_cat,alpha_ahp,lambda,n_oxyr,k_oxyr_red,k_oxyr_ox,Km_oxyr);
f_as = x_as(:,2).^n_oxyr./(x_as(:,2).^n_oxyr*(1+k_oxyr_red/k_oxyr_ox)+k_oxyr_red/k_oxyr_ox*Km_oxyr^n_oxyr);

%% plot Hout vs Fold change
%figure();
xbound = [-0.01,0.5];

subplot(5,1,1);
hold on;
plot([t_bs-max(t_bs);t_as]/3600,[x_bs(:,1);x_as(:,1)]*1e6);
xlim(xbound);
box on;

subplot(5,1,2);
hold on;
plot([t_bs-max(t_bs);t_as]/3600,[x_bs(:,2);x_as(:,2)]*1e6);
xlim(xbound);
box on;

subplot(5,1,3);
hold on;
plot([t_bs-max(t_bs);t_as]/3600,[f_bs;f_as]);
xlim(xbound);
box on;

subplot(5,1,4);
hold on;
plot([t_bs-max(t_bs);t_as]/3600,[x_bs(:,3);x_as(:,3)]*1e6);
xlim(xbound);
box on;

subplot(5,1,5);
hold on;
plot([t_bs-max(t_bs);t_as]/3600,[x_bs(:,4);x_as(:,4)]*1e6);
xlim(xbound);
box on;
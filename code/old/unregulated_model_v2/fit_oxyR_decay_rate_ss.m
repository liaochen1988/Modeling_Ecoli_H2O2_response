clear all;
clc;

%% oxyR oxidation rate at different h2o2 concentration
xdata_h2o2_kobs = [2.976,5.107,7.167,9.582,11.982,14.057,15.945,19.499,21.934,...
            27.088,32.786,44.054,65.975,88.398,110.454,132.145,154.433,176.341,198.563]; % uM

ydata_h2o2_kobs = [0.300,0.679,0.989,1.265,1.556,1.887,2.246,2.654,2.871,3.416,3.899,5.451,6.644,...
         7.232,8.178,8.233,8.708,8.865,8.595]; % 1/s
     
%% oxidized oxyR fraction at different h2o2 concentration
xdata_h2o2_foxyrox = [0,0.025,0.05,0.075,0.1,0.2,0.5,1.0,2.0,5.0,10.0]; % uM
ydata_h2o2_foxyrox = [0,0.216,0.300,0.467,0.533,0.732,0.865,0.937,0.966,0.978,1.000]; % percentage

%% oxidized oxyR fraction at different GSH/GSSG ratio
xdata_gshgssgr_foxyrox = [1.944e-4,0.001,0.005,0.011,0.020,0.031,0.044,0.124];
ydata_gshgssgr_foxyrox = [99.925,93.740,92.639,40.366,5.558,3.062,1.083,0.708]/100;

%% fit a Hill function to oxyR oxidation rate first
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [1,50,10]; % Km,n,k_obs
lb = [0,0,0];
ub = [Inf,Inf,Inf];
[sol,~,~,exitflag] = lsqnonlin(@mse_vary_h2o2_kobs,p0,lb,ub,options,xdata_h2o2_kobs,ydata_h2o2_kobs);
assert(exitflag>0);

n_oxyr = sol(1);
Km_oxyr = sol(2);
k_oxyr_ox = sol(3);

%% fit a Hill function to oxyR fraction
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [0]; % Keq
lb = [-Inf];
ub = [Inf];
[sol,~,~,exitflag] = lsqnonlin(@mse_vary_gshgssgr_foxyrox,p0,lb,ub,options,xdata_gshgssgr_foxyrox,ydata_gshgssgr_foxyrox);
assert(exitflag>0);

Keq = 10^sol(1);

%% fit a Hill function to oxyR fraction
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [-4]; % k_oxyr_ox1
lb = [-Inf];
ub = [Inf];
[sol,~,~,exitflag] = lsqnonlin(@mse_vary_h2o2_foxyrox_ss,p0,lb,ub,options,xdata_h2o2_foxyrox,ydata_h2o2_foxyrox,n_oxyr,Km_oxyr,k_oxyr_ox,Keq);
assert(exitflag>0);

k_oxyr_ox1 = 10^sol(1);
k_oxyr_red1 = Keq*k_oxyr_ox1;

%% plot fitting
figure();

subplot(1,3,1);
hold on;

xfit = 10.^[0:0.1:3];
yfit = k_oxyr_ox*xfit.^n_oxyr./(xfit.^n_oxyr+Km_oxyr^n_oxyr);
plot(xdata_h2o2_kobs,ydata_h2o2_kobs,'ko');
plot(xfit,yfit,'k-');
set(gca,'XScale','log');
axis square;
box on;
axis([1e0,1e3,0,10]);
xlabel('[H_2O_2] (\muM)');
ylabel('H2O2-induced OxyR oxdization rate (s^{-1})');

subplot(1,3,2);
hold on;
xfit = 10.^[-4:0.1:0]; % GSH^2/GSSG
yfit = 1./(1+Keq*xfit.^4);

plot(xdata_gshgssgr_foxyrox,ydata_gshgssgr_foxyrox,'ko');
plot(xfit,yfit,'k-');
set(gca,'XScale','log');
axis square;
box on;
xlabel('[GSH]^2/[GSSG] (M)');
ylabel('Oxidized OxyR fraction');
axis([1e-4,1e0,-0.1,1.1]);

subplot(1,3,3);
hold on;
xfit = 10.^[-4:0.1:2];

GSSG = 1e-4; % M
GSH = 0.025; % M
kswitch = k_oxyr_ox*xfit.^n_oxyr./(Km_oxyr^n_oxyr+xfit.^n_oxyr);
yfit = (kswitch+k_oxyr_ox1*GSSG.^4)./(k_oxyr_red1*GSH.^8+kswitch+k_oxyr_ox1*GSSG.^4);

plot(xdata_h2o2_foxyrox,ydata_h2o2_foxyrox,'ko');
plot(xfit,yfit,'k-');
set(gca,'XScale','log');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');
ylabel('Oxidized OxyR fraction');
axis([1e-4,1e2,-0.1,1.1]);
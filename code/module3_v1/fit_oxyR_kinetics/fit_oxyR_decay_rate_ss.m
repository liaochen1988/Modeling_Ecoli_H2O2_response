clear all;
clc;

%% fit a Hill function to oxyR oxidation rate
% Lee2004NatStrucMolBiol

h2o2_lee2004 = [2.976,5.107,7.167,9.582,11.982,14.057,15.945,19.499,21.934,...
            27.088,32.786,44.054,65.975,88.398,110.454,132.145,154.433,176.341,198.563]; % uM

kobs_lee2004 = [0.300,0.679,0.989,1.265,1.556,1.887,2.246,2.654,2.871,3.416,3.899,5.451,6.644,...
         7.232,8.178,8.233,8.708,8.865,8.595]; % 1/s
     
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [1,50,10]; % n_oxyr, Km_oxyr, kon_oxyr_h2o2
lb = [0,0,0];
ub = [Inf,Inf,Inf];
[sol,~,~,exitflag] = lsqnonlin(@mse_vary_h2o2_kobs,p0,lb,ub,options,h2o2_lee2004,kobs_lee2004);
assert(exitflag>0);

n_oxyr = sol(1);
Km_oxyr = sol(2);
kon_oxyr_h2o2 = sol(3);

%% fit a Hill function to oxidized oxyR fraction at different GSH^2/GSSG ratio
% Zheng1998Science

gsh2gssg_zheng1998 = [1.944e-4,0.001,0.005,0.011,0.020,0.031,0.044,0.124];
frac_oxyr_ox_zheng1998 = [99.925,93.740,92.639,40.366,5.558,3.062,1.083,0.708]/100;

options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [0]; % Keq
lb = [-Inf];
ub = [Inf];
[sol,~,~,exitflag] = lsqnonlin(@mse_vary_gshgssgr_foxyrox,p0,lb,ub,options,gsh2gssg_zheng1998,frac_oxyr_ox_zheng1998);
assert(exitflag>0);

Keq = 10^sol(1); % koff_oxyr_grx/kon_oxyr_grx

%% fit a Hill function to oxyR fraction
% Aslund1999PNAS

h2o2_aslund1999 = [0,0.025,0.05,0.075,0.1,0.2,0.5,1.0,2.0,5.0,10.0]; % uM
frac_oxyr_ox_aslund1999 = [0,0.216,0.300,0.467,0.533,0.732,0.865,0.937,0.966,0.978,1.000]; % percentage

options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [0]; % kon_oxyr_grx
[sol,~,~,exitflag] = lsqnonlin(@mse_vary_h2o2_foxyrox_ss,p0,[],[],options,h2o2_aslund1999,frac_oxyr_ox_aslund1999,n_oxyr,Km_oxyr,kon_oxyr_h2o2,Keq);
assert(exitflag>0);

kon_oxyr_grx = 10^sol(1);
koff_oxyr_grx = Keq*kon_oxyr_grx;

%% plot fitting
figure();

subplot(1,3,1);
hold on;

xfit = 10.^[0:0.1:3];
yfit = kon_oxyr_h2o2*xfit.^n_oxyr./(xfit.^n_oxyr+Km_oxyr^n_oxyr);
plot(h2o2_lee2004,kobs_lee2004,'ko');
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

plot(gsh2gssg_zheng1998,frac_oxyr_ox_zheng1998,'ko');
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
Grx = 10*1e-6; % M
kswitch = kon_oxyr_h2o2*xfit.^n_oxyr./(Km_oxyr^n_oxyr+xfit.^n_oxyr);
yfit = (kswitch+kon_oxyr_grx*Grx*GSSG.^4)./(koff_oxyr_grx*Grx*GSH.^8+kswitch+kon_oxyr_grx*Grx*GSSG.^4);

plot(h2o2_aslund1999,frac_oxyr_ox_aslund1999,'ko');
plot(xfit,yfit,'k-');
set(gca,'XScale','log');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');
ylabel('Oxidized OxyR fraction');
axis([1e-4,1e2,-0.1,1.1]);
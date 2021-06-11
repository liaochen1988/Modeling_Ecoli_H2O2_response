clear all;
clc;

%% final values
% hill_c = 1.3614
% Km = 41.3179
% kobs_max= 9.9863 1/s
% kd = 0.0023 1/s, corresponding to half life 5 min

%% oxyR oxidation rate vs H2O2
h2o2_swtich = [2.976,5.107,7.167,9.582,11.982,14.057,15.945,19.499,21.934,...
            27.088,32.786,44.054,65.975,88.398,110.454,132.145,154.433,176.341,198.563]; % uM

switch_rate = [0.300,0.679,0.989,1.265,1.556,1.887,2.246,2.654,2.871,3.416,3.899,5.451,6.644,...
         7.232,8.178,8.233,8.708,8.865,8.595]; % 1/s
     
%% oxidized oxyR fraction vs H2O2
h2o2_foxyrox = [0,0.025,0.05,0.075,0.1,0.2,0.5,1.0,2.0,5.0,10.0]; % uM
oxyr_ox_frac1 = [0,0.216,0.300,0.467,0.533,0.732,0.865,0.937,0.966,0.978,1.000]; % percentage

% plot(h2O2_obs,k_obs,h2o2_invitro*1e6,oxyr_ox_invitro);

%% oxidized oxyR fraction vs GSH/GSSG ratio
GSHGSSGR = [1.944e-4,0.001,0.005,0.011,0.020,0.031,0.044,0.124];
oxyr_ox_frac2 = [99.925,93.740,92.639,40.366,5.558,3.062,1.083,0.708]/100;

%% fit a Hill function to oxyR oxidation rate first
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [1,50,10]; % Km,n,k_obs
lb = [0,0,0];
ub = [Inf,Inf,Inf];
[hill_opt_1,~,~,exitflag] = lsqnonlin(@mseHill_v2,p0,lb,ub,options,h2o2_swtich,switch_rate);
assert(exitflag>0);

n_oxyr = hill_opt_1(1);
Km_oxyr = hill_opt_1(2);
k_oxyr_ox = hill_opt_1(3);

%% fit a Hill function to oxyR fraction
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [0]; % Keq
lb = [-Inf];
ub = [Inf];
[hill_opt_2,~,~,exitflag] = lsqnonlin(@mseHill_v5,p0,lb,ub,options,GSHGSSGR,oxyr_ox_frac2);
assert(exitflag>0);

Keq = 10^hill_opt_2(1);

%% fit a Hill function to oxyR fraction
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [-20]; % k_oxyr_ox1
lb = [-Inf];
ub = [Inf];
[hill_opt_3,~,~,exitflag] = lsqnonlin(@mseHill_v4,p0,lb,ub,options,h2o2_foxyrox,oxyr_ox_frac1,hill_opt_1,Keq);
assert(exitflag>0);

k_oxyr_ox1 = 10^hill_opt_3(1);
k_oxyr_red1 = Keq*k_oxyr_ox1;

%% plot fitting
figure();

%subplot(1,3,1);
hold on;

xfit = 10.^[-4:0.1:4];
yfit = k_oxyr_ox*xfit.^n_oxyr./(xfit.^n_oxyr+Km_oxyr^n_oxyr);
yyaxis left
plot(h2o2_swtich,switch_rate,'ko');
plot(xfit,yfit,'k-');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');
ylabel('k_{obs} (s^{-1})');

%subplot(1,3,2);
hold on;
xfit = 10.^[-4:0.1:4]; % GSH^2/GSSG
yfit = 1./(1+Keq*xfit.^4);

yyaxis right
plot(GSHGSSGR,oxyr_ox_frac2,'ko');
plot(xfit,yfit,'k-');
set(gca,'XScale','log');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');

%subplot(1,3,3);
hold on;
xfit = 10.^[-4:0.1:4];

GSSG = 0.1; % mM
GSH = 25; % mM
kswitch = k_oxyr_ox*xfit.^n_oxyr./(Km_oxyr^n_oxyr+xfit.^n_oxyr);
yfit = (kswitch+k_oxyr_ox1*GSSG.^4)./(k_oxyr_red1*GSH.^8+kswitch+k_oxyr_ox1*GSSG.^4);

yyaxis right
plot(h2o2_foxyrox,oxyr_ox_frac1,'ko');
plot(xfit,yfit,'k-');
set(gca,'XScale','log');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');
ylabel('OxyR-ox fraction');
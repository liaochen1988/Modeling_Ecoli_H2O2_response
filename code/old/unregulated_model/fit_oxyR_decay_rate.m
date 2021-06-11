clear all;
clc;

%% final values
% hill_c = 1.3614
% Km = 41.3179
% kobs_max= 9.9863 1/s
% kd = 0.0023 1/s, corresponding to half life 5 min

%% oxyR oxidation rate
h2o2_invitro_1 = [2.976,5.107,7.167,9.582,11.982,14.057,15.945,19.499,21.934,...
            27.088,32.786,44.054,65.975,88.398,110.454,132.145,154.433,176.341,198.563]; % uM

oxidation_rate = [0.300,0.679,0.989,1.265,1.556,1.887,2.246,2.654,2.871,3.416,3.899,5.451,6.644,...
         7.232,8.178,8.233,8.708,8.865,8.595]; % 1/s
     
%% oxidized oxyR fraction
h2o2_invitro_2 = [0,0.025,0.05,0.075,0.1,0.2,0.5,1.0,2.0,5.0,10.0]; % uM
oxidation_fraction = [0,0.216,0.300,0.467,0.533,0.732,0.865,0.937,0.966,0.978,1.000]; % percentage

% plot(h2O2_obs,k_obs,h2o2_invitro*1e6,oxyr_ox_invitro);

%% fit a Hill function to oxyR oxidation rate first
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [1,50,10]; % Km,n,k_obs
lb = [0,0,0];
ub = [Inf,Inf,Inf];
[hill_opt_1,~,~,exitflag] = lsqnonlin(@mseHill_v2,p0,lb,ub,options,h2o2_invitro_1,oxidation_rate);
assert(exitflag>0);

hill_c = hill_opt_1(1);
Km = hill_opt_1(2);
kobs_max = hill_opt_1(3);

%% fit a Hill function to oxyR fraction
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [1e-3]; % Km,n,k_obs
lb = [0];
ub = [Inf];
[hill_opt_2,~,~,exitflag] = lsqnonlin(@mseHill_v3,p0,lb,ub,options,h2o2_invitro_2,oxidation_fraction,hill_opt_1);
assert(exitflag>0);

kd = hill_opt_2(1);

%% plot fitting
figure();

%subplot(1,2,1);
hold on;

xfit = 10.^[-2:0.1:3];
yfit = kobs_max*xfit.^hill_c./(xfit.^hill_c+Km^hill_c);
yyaxis left
plot(h2o2_invitro_1,oxidation_rate,'ko');
plot(xfit,yfit,'k-');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');
ylabel('k_{obs} (s^{-1})');

%subplot(1,2,2);
hold on;

xfit = 10.^[-2:0.1:3];
yfit = xfit.^hill_c./(xfit.^hill_c*(1+kd/kobs_max)+kd/kobs_max*Km^hill_c);
yyaxis right
plot(h2o2_invitro_2,oxidation_fraction,'ko');
plot(xfit,yfit,'k-');
set(gca,'XScale','log');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');
ylabel('OxyR-ox fraction');
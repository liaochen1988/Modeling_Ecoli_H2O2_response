clear all;
clc;

%% parameters
P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 0;
Km_ahp = 1.2e-6; % M
k_cat = 2.7e-13; % L/s
k_ahp = 2.1e-18; % mol/s

%% oxyR oxidation data
h2o2_invitro = [0,0.025,0.05,0.075,0.1,0.2,0.5,1.0,2.0,5.0,10.0]/1e6; % M
oxyr_ox_invitro = [0,0.216,0.300,0.467,0.533,0.732,0.865,0.937,0.966,0.978,1.000];

h2o2_invivo = [0,0.5,1.0,2.0,5.0,10.0]/1e6; % M
oxyr_ox_invivo = [0,0,0,0,0.6,1.0];

%% fit Hill function
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [1.5,0.1];
lb = [0,0];
ub = [Inf,Inf];
[p_opt_invitro,~,~,exitflag] = lsqnonlin(@mseHill,p0,lb,ub,options,h2o2_invitro,oxyr_ox_invitro);
assert(exitflag>0);

p0 = [10.7,4.8e-6];
lb = [0,0];
ub = [Inf,Inf];
[p_opt_invivo,~,~,exitflag] = lsqnonlin(@mseHill,p0,lb,ub,options,h2o2_invivo,oxyr_ox_invivo);
assert(exitflag>0);

%% Predict OxyR oxidation percentage
Hout = 10.^[-9:0.1:-4]; % M
Delta = Km_ahp*PA+Km_ahp*k_cat-Hout*PA-k_met+k_ahp;
Hin = (-Delta+sqrt(Delta.^2+4*Km_ahp*(PA+k_cat)*(k_met+Hout*PA)))/2/(PA+k_cat);
OxyR = Hin.^p_opt_invitro(1)./(Hin.^p_opt_invitro(1)+p_opt_invitro(2)^p_opt_invitro(1));

%% plot Hout vs OxyR
figure();

subplot(1,2,1);
hold on;
plot(Hout*1e6, OxyR, 'r-');
plot(Hin*1e6, OxyR,'b-');
legend('In vivo','In vitro');
title('Simulation');
axis([1e-3,1e2,-0.1,1.1]);
axis square;
box on;
ylabel('Oxidized OxyR (%)');
xlabel('Extracellular H_2O_2 (\muM)');
set(gca,'XScale','log');
set(gca,'XTick',10.^[-3:1:2]);
set(gca,'YTick',[0:0.2:1.0]);

subplot(1,2,2);
hold on;

xfit = 10.^[-3:0.01:2]/1e6;
yfit_invitro = xfit.^p_opt_invitro(1)./(xfit.^p_opt_invitro(1)+p_opt_invitro(2)^p_opt_invitro(1));
yfit_invivo = xfit.^p_opt_invivo(1)./(xfit.^p_opt_invivo(1)+p_opt_invivo(2)^p_opt_invivo(1));

plot(h2o2_invivo*1e6, oxyr_ox_invivo, 'ko','MarkerFaceColor','r');
plot(xfit*1e6,yfit_invivo,'r-');
plot(h2o2_invitro*1e6, oxyr_ox_invitro, 'ko','MarkerFaceColor','b');
plot(xfit*1e6,yfit_invitro,'b-');
legend('In vivo','In vitro');
title('Measurement');
axis([1e-3,1e2,-0.1,1.1]);
axis square;
box on;
ylabel('Oxidized OxyR (%)');
xlabel('Extracellular H_2O_2 (\muM)');
set(gca,'XScale','log');
set(gca,'XTick',10.^[-3:1:2]);
set(gca,'YTick',[0:0.2:1.0]);



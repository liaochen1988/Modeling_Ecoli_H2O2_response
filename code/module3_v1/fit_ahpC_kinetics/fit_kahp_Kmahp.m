clear all;
clc;

%% Experimental data used in parameter fitting

% (In vitro) Change of oxidized OxyR fraction at varied H2O2 concentration
% Aslund1999PNAS
h2o2_aslund1999_invitro = [0,0.025,0.05,0.075,0.1,0.2,0.5,1.0,2.0,5.0,10.0]/1e6; % M
frac_oxyr_ox_aslund1999_invitro = [0,0.216,0.300,0.467,0.533,0.732,0.865,0.937,0.966,0.978,1.000];

% (In vivo) Change of oxidized OxyR fraction at varied H2O2 concentration
% Aslund1999PNAS
h2o2_aslund1999_invivo = [0,0.5,1.0,2.0,5.0,10.0]/1e6; % M
frac_oxyr_ox_aslund1999_invivo = [0,0,0,0,0.6,1.0];

% (In vivo) H2O2 decomposition rate at varied H2O2 concentration
% SeaverImlay2001JB
h2o2_seaver2001 = [0,0.585,0.875,1.368,4.942,4.967,5.007,9.961,15.098,15.023,27.982,39.990,39.976]/1e6; %M
Jahp_seaver2001 = [0.024,1.398,2.176,3.023,4.085,5.055,9.491,13.882,14.477,16.730,15.500,16.667,17.704]; % uM/min*OD

%% Fit k_ahp and Km_ahp
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-12,'TolFun',1e-12,'MaxFunEval',Inf,'Algorithm','levenberg-marquardt');
p0 = [log10(2.1e-18),log10(1.2e-6)]; % use values from SeaverImlay2001JB for initial guess 
lb = [-Inf,-Inf];
ub = [0,0];
[sol,~,~,exitflag] = lsqnonlin(@mse_vary_h2o2,p0,lb,ub,options,h2o2_aslund1999_invivo,frac_oxyr_ox_aslund1999_invivo,h2o2_seaver2001,Jahp_seaver2001);
assert(exitflag>0);

k_ahp = 10.^(sol(1));
Km_ahp = 10.^(sol(2));

%% The fitted value k_ahp = 2.95e-18, which is similar to the reported value 2.1e-18 in Seaver2001JB

%% Plot comparison between experiment and simulation
[~,Hout,Hin_f,frac_oxyr_ox_pred_invivo,Hin_J,Jahp_pred] = ...
    mse_vary_h2o2(sol,h2o2_aslund1999_invivo,frac_oxyr_ox_aslund1999_invivo,h2o2_seaver2001,Jahp_seaver2001);

%figure();

subplot(1,2,1);
hold on;

% calculate in vitro prediction of oxidized oxyR fraction
load('par_est_oxyr_kinetics');
Hout_invitro = 10.^[-4:0.1:2]; % M
GSSG = 1e-4; % M
GSH = 0.025; % M
Grx = 10*1e-6; % M
kswitch = kon_oxyr_h2o2*Hout_invitro.^n_oxyr./(Km_oxyr^n_oxyr+Hout_invitro.^n_oxyr);
frac_oxyr_ox_pred_invitro = (kswitch+kon_oxyr_grx*Grx*GSSG.^4)./(kswitch+koff_oxyr_grx*Grx*GSH.^8+kon_oxyr_grx*Grx*GSSG.^4);

plot(h2o2_aslund1999_invitro*1e6, frac_oxyr_ox_aslund1999_invitro, 'ko','MarkerFaceColor','b');
plot(Hout_invitro, frac_oxyr_ox_pred_invitro,'b-');
plot(h2o2_aslund1999_invivo*1e6, frac_oxyr_ox_aslund1999_invivo, 'ko','MarkerFaceColor','r');
plot(Hout*1e6, frac_oxyr_ox_pred_invivo, 'r-');
legend('In vitro obs','In vitro sim','In vivo obs','In vivo sim');
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
plot(Hout*1e6,Jahp_pred/max(Jahp_pred),'k-');
plot(h2o2_seaver2001*1e6,Jahp_seaver2001/max(Jahp_seaver2001),'ko');
axis square;
box on;

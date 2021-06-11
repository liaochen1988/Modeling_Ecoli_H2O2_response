clear all;
clc;

%% Fixed parameters

% sugar oxidation rate and Michaelis constant 
k_o = 0.0508; % mol/L/h
Km_o = 0.0168; % mol/L

% H2O2 removal rate and Michaelis constant
k_d = 1.7884e+04; % umol/L/h
Km_d = 5.9e-3*1e6; % umol/L

% H2O2 production rate per carbon uptake rate
f_u = 6.4799e3;

% H2O2 production rate per carbon oxidation rate
f_o = 35.7342; % convert from carbon consumed (M) to h2o2 produced (uM), i.e., metabolizing 1M C produces 5.165 uM H2O2

% Vin
Vin = 1e-15; % L

% Vout
Vout = 0.2e-3; % L

% switch time
phi = 0.3085; % hour

% sugar uptake rate and Michaelis constant

% Transport of Glycerol by Pseudomonas aeruginosa,SAN-SAN TSAY,1971
% 21.3 nmol/min/(mg dry wt cell) = 151e-6*60*1000 = 9.06 mmol/gdw/h

% The effect of nutrient limitation on glycerol uptake and metabolism in continuous cultures of Pseudornonas aeruginosa
% 21.3 nmol/min/(mg dry wt cell) = 21.3e-6*60*1000 = 1.3 mmol/gdw/h

% glycerol: C3H8O3
k_u = 4.2114; % mol/L/h
Km_u = 4.8e-4; % mol/L

% half-maixmal inhibition constant and Hill coefficient of H2O2
Ki_u = 0.5; % uM, close to minimum value which fits equally well
n_u = 1;

% bacterial growth rate per carbon uptake rate
% PRODUCTION OF BIOSURFACTANTS FROM Pseudomonas aeruginosa PA1 ISOLATED IN OIL ENVIRONMENTS
% complete consumption of 15 g/L glycerol results in 2.5 g/L biomass
e_g = 0.0033;

% diffusion rate
% A Kinetic Platform to Determine the Fate of Hydrogen Peroxide in Escherichia coli
% it is estimated to be the lower bound
k_diff = 0.25*3600; % 1/h

% basal h2o2 production rate
r_g0 = 0.025;

% maintenance cost
r_mc = 0.015;

%% Cell growth model parameters
beta = 3.0e6; % uM
m_r = 11738;
m_e = 325;
k_r = 7.56e4; % 1/h
K_ma = 20; % uM
K_ix = 60; % uM
k_x = 3.6e3; % 1/h
d_x = 1.26e2; % 1/hr
K_ia = 1e4; % uM

%% Parameter fitting
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);

rawDataEM = xlsread('./181130_glucose_Amplex.xlsx','AmplexEM');
timeEM = rawDataEM(:,2)/3600;  % in unit of hour
rawDataEM = rawDataEM(:,[4:end]);
h2o2EM = mean(rawDataEM(:,[10:12]),2)/1189.6; % uM

rawDataOD = xlsread('./181130_glucose_Amplex.xlsx','OD600');
timeOD = rawDataOD(:,2)/3600;  % in unit of hour
rawDataOD = rawDataOD(:,[4:end]);
h2o2OD = mean(rawDataOD(:,[10:12])-rawDataOD(:,[2:4]),2)*5e8;

initE = 1000;
p0 = [k_d,k_u,e_g,phi,initE];
lb = [0,0,0,0,0];
ub = [Inf,Inf,Inf,1,beta/m_e];
[p_opt,~,~,exitflag] = lsqnonlin(@mse_h2o2_removal_model_R,p0,lb,ub,options,timeOD,h2o2OD,timeEM,h2o2EM,...
    k_o,Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,k_diff,Vin,Vout,r_g0,r_mc,beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia);
assert(exitflag>0);

k_d = p_opt(1);
k_u = p_opt(2);
e_g = p_opt(3);
phi = p_opt(4);
initE = p_opt(5);

%% Simulation

% k_d > 0
x0 = [0, 3.44, 3/12, 1.1e7, 0, 1, 0, initE];
tspan = [0,48];
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,8),'NonNegative',1:8);
[t_wc,x_wc] = ode15s(@h2o2_removal_model_R,tspan,x0,options,...
    phi,k_o,k_u,k_d,Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,e_g,k_diff,Vin,Vout,r_g0,r_mc,...
    beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia);

% no cell
x0(2) = 4.1;
x0(4)=0;
[t_nc,x_nc] = ode15s(@h2o2_removal_model_R,tspan,x0,options,...
    phi,k_o,k_u,k_d,Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,e_g,k_diff,Vin,Vout,r_g0,r_mc,...
    beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia);

%% plot
figure();

% Intracellular H2O2
subplot(2,2,1);
hold on;
plot(t_wc,x_wc(:,1),'r-');
plot(t_nc,x_nc(:,1),'b-');
xlabel('time (h)');
ylabel('In h2o2 (\muM)');
xlim([0,20]);

rawDataEM = xlsread('./181130_glucose_Amplex.xlsx','AmplexEM');
time_obs = rawDataEM(:,2)/3600;  % in unit of hour
rawDataEM = rawDataEM(:,[4:end]);
h2o2_obs_nc = mean(rawDataEM(:,[2:4]),2)/1189.6; % uM
h2o2_obs_wc = mean(rawDataEM(:,[10:12]),2)/1189.6; % uM

% Extracellular H2O2
subplot(2,2,2);
hold on;
plot(t_wc,x_wc(:,2),'r-');
plot(t_nc,x_nc(:,2),'b-');
plot(time_obs,h2o2_obs_wc,'r--');
plot(time_obs,h2o2_obs_nc,'b--');
xlabel('time (h)');
ylabel('Ex h2o2 (\muM)');
xlim([0,20]);

% Carbon

subplot(2,2,3);
hold on;
plot(t_wc,x_wc(:,3),'r-');
plot(t_nc,x_nc(:,3),'b-');
xlabel('time (h)');
ylabel('carbon (mM)');
xlim([0,20]);

% Cell number

rawDataOD = xlsread('./181130_glucose_Amplex.xlsx','OD600');
time_obs = rawDataOD(:,2)/3600;  % in unit of hour
rawDataOD = rawDataOD(:,[4:end]);
celln_obs_wc = mean(rawDataOD(:,[10:12])-rawDataOD(:,[2:4]),2)*5e8;

subplot(2,2,4);
hold on;
plot(t_wc,x_wc(:,4),'r-');
plot(t_nc,x_nc(:,4),'b-');
plot(time_obs,celln_obs_wc,'r--');
xlabel('time (h)');
ylabel('cell number');
xlim([0,20]);
set(gca,'YScale','log');
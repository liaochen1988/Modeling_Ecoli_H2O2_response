clear all;
clc;

addpath('../apm');

%% Dynamic optimization

% Select server
server = 'http://byu.apmonitor.com';
%server = 'http://xps.apmonitor.com';

%   objective: total biomass production
app = 'cellgrowth_total_biomass';

%   write data file
%write_data_file_total_biomass([0:0.1:0.7,0.8:0.2:2,2.5:0.5:3.5,3.6:0.2:5.8,6:0.5:10]);
write_data_file_total_biomass([0,10.^[-6:-3],0.01:0.001:0.5,0.51:0.01:4]);

% Clear previous application
apm(server,app,'clear all');

% Load model
apm_load(server,app,[app '.apm']);

% Load csv file
csv_load(server,app,[app '.csv']);

% Option to select
apm_option(server,app,'nlc.nodes',2);
apm_option(server,app,'nlc.solver',3);
apm_option(server,app,'nlc.imode',6);

apm_info(server,app,'MV','phi');
apm_option(server,app,'phi.status',1);
apm_option(server,app,'apm.max_iter',10000);
apm_option(server,app,'phi.dcost',0);
apm_option(server,app,'apm.otol',1e-7);
apm_option(server,app,'apm.rtol',1e-7);
apm_option(server,app,'phi.lower',0);
apm_option(server,app,'phi.upper',1);

% Solve on APM server
output = apm(server,app,'solve');
disp(output);
success_check = apm_tag(server,'cellgrowth_total_biomass','apm.appstatus');
if(~success_check)
    apm_get(server,app,'infeasibilities.txt');
    error('Solution not found by APMonitor!');
end

% Retrieve results
y = apm_sol(server,app);
optimalSol = y.x;

%% plot phi(t)

figure();
hold on;

plot(optimalSol.time,optimalSol.phi,'-k');
 
axis square;
box on;
xlabel('Time (h)');
ylabel('phi');

save('optimalSol');
 
%% plot simulation

load('optimalSol');

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
r_mc = 0; % 0.015;

beta = 3.0e6; % uM
m_r = 11738;
m_e = 325;
k_r = 7.56e4; % 1/h
K_ma = 20; % uM
K_ix = 60; % uM
k_x = 3.6e3; % 1/h
d_x = 1.26e2; % 1/hr
K_ia = 1e4; % uM

% k_d > 0
x0 = [0, 3.44, 3/12, 1.1e7, 1, 1, 1, 1000];
tspan = [0,5];
options = odeset('RelTol',1e-6,'AbsTol',1e-6*ones(1,8),'NonNegative',1:8);
[t_wt,x_wt] = ode15s(@h2o2_removal_model_R,tspan,x0,options,...
    phi,k_o,k_u,k_d,Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,e_g,k_diff,Vin,Vout,r_g0,r_mc,...
    beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia);
[t_opt_dyn,x_opt_dyn] = ode15s(@h2o2_removal_model_Rs,tspan,x0,options,...
    optimalSol,k_o,k_u,k_d,Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,e_g,k_diff,Vin,Vout,r_g0,r_mc,...
    beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia);    

figure();

% Intracellular H2O2
subplot(2,3,1);
hold on;
plot(t_wt,x_wt(:,1),'r-');
plot(t_opt_dyn,x_opt_dyn(:,1),'b-');
xlabel('time (h)');
ylabel('In h2o2 (\muM)');
axis square;
box on;
%axis([0,10,0,1.5]);

% Extracellular H2O2

rawDataEM = xlsread('./181130_glucose_Amplex.xlsx','AmplexEM');
time_obs = rawDataEM(:,2)/3600;  % in unit of hour
rawDataEM = rawDataEM(:,[4:end]);
h2o2_obs_nc = mean(rawDataEM(:,[2:4]),2)/1189.6; % uM
h2o2_obs_wc = mean(rawDataEM(:,[10:12]),2)/1189.6; % uM

subplot(2,3,2);
hold on;
plot(t_wt,x_wt(:,2),'r-');
plot(time_obs,h2o2_obs_wc,'r.');
plot(t_opt_dyn,x_opt_dyn(:,2),'b-');
xlabel('time (h)');
ylabel('Ex h2o2 (\muM)');
axis square;
box on;
%axis([0,10,0,1.2]);

% Carbon

subplot(2,3,3);
hold on;
plot(t_wt,x_wt(:,3),'r-');
plot(t_opt_dyn,x_opt_dyn(:,3),'b-');
xlabel('time (h)');
ylabel('carbon (mM)');
axis square;
box on;
%axis([0,10,0,0.25]);

% Cell number

rawDataOD = xlsread('./181130_glucose_Amplex.xlsx','OD600');
time_obs = rawDataOD(:,2)/3600;  % in unit of hour
rawDataOD = rawDataOD(:,[4:end]);
celln_obs_wc = mean(rawDataOD(:,[10:12])-rawDataOD(:,[2:4]),2)*5e8;

subplot(2,3,4);
hold on;
plot(t_wt,x_wt(:,4),'r-');
plot(time_obs,celln_obs_wc,'r.');
plot(t_opt_dyn,x_opt_dyn(:,4),'b-');
xlabel('time (h)');
ylabel('cell number');
set(gca,'YScale','log');
axis square;
box on;
axis([0,10,1e7,1e9]);

subplot(2,3,5);
hold on;
plot(t_wt,ones(1,length(t_wt))*phi,'r-');
plot(optimalSol.time,optimalSol.phi,'b-');
xlabel('time (h)');
ylabel('control function');
axis square;
box on;
axis([0,10,0,1]);

%% plot

figure();

subplot(1,3,1);
plot(optimalSol.time,optimalSol.phi,'-k');
axis square;
box on;
xlabel('Time (h)');
ylabel('phi');
axis([0,4,-0.1,1.1]);

subplot(1,3,2);
hold on;
plot(t_opt_dyn,x_opt_dyn(:,6)*m_r/beta,'r-');
plot(t_opt_dyn,x_opt_dyn(:,8)*m_e.*pchip(optimalSol.time,optimalSol.phi/beta,t_opt_dyn),'b-');
plot(t_opt_dyn,x_opt_dyn(:,8)*m_e.*(1-pchip(optimalSol.time,optimalSol.phi/beta,t_opt_dyn))/beta,'g-');
axis square;
box on;
xlabel('time (h)');
ylabel('proteome fraction');
axis([0,4,-0.1,1.1]);

subplot(1,3,3);
hold on;
plot(t_opt_dyn,x_opt_dyn(:,4),'k-');
xlabel('time (h)');
ylabel('cell number');
set(gca,'YScale','log');
axis square;
box on;
axis([0,4,1e7,1e9]);

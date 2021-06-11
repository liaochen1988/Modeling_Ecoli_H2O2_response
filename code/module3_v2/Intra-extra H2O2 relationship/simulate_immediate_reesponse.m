clear all;
clc;

P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20; % mol/s
Vc = 3.23*1e-15; % L
Ve = 0.2e-3;  % L

k_kat = 1.63e4; % 1/s
Km_kat = 3.9/1e3;% M
Ki_h2o2 = 19.55/1e6; %M
alpha_max_kat = 19.8/1e6/3600; % M/s

k_ahp = 52.4; % 1/s
Km_ahp = 1.4987e-07; % M
alpha_max_ahp = 12.54/1e6/3600; % M/s

alpha_max_grx = 19.39/1e6/3600; % M/s

n_oxyr = 1.3614;
Km_oxyr = 41.3179e-6; % M
kon_oxyr_h2o2 = 9.99; % 1/s
kon_oxyr_grx = 1.49e7; % 1/s/M^5
koff_oxyr_grx = 1.52e15; % 1/s/M^9

GSSG = 1e-4; % M
GSH = 0.025; % M
Grx = 10*1e-6; % M
lambda = 0.99/3600; % 1/s

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,5),'NonNegative',1:5);

% pre-addition
x0 = [0,1,1,1,1];
N = 0;
[~,y0] = ode15s(@h2o2_model,[0,1e10],x0,options,N,GSH,GSSG,lambda,...
                   PA,k_met,Vc,Ve,...
                   k_kat,Km_kat,Ki_h2o2,alpha_max_kat,...
                   k_ahp,Km_ahp,alpha_max_ahp,...
                   alpha_max_grx,...
                   n_oxyr,Km_oxyr,kon_oxyr_h2o2,kon_oxyr_grx,koff_oxyr_grx);
    
% after addition
yic = y0(end,:);
yic(1) = 10 * 1e-6; % M
yic(3) = 0;
yic(4) = 0;
N = 0.01 * 5e8;
[t1,y1] = ode15s(@h2o2_model,[0,3600],yic,options,N,GSH,GSSG,lambda,...
                   PA,k_met,Vc,Ve,...
                   k_kat,Km_kat,Ki_h2o2,alpha_max_kat,...
                   k_ahp,Km_ahp,alpha_max_ahp,...
                   alpha_max_grx,...
                   n_oxyr,Km_oxyr,kon_oxyr_h2o2,kon_oxyr_grx,koff_oxyr_grx);
               
%% plot

figure();

subplot(2,2,1);
hold on;
plot(t1/3600, y1(:,1)*1e6);
axis([0,1,0,10]);
xlabel('Time (h)');
ylabel('[H_o] (\muM)');
axis square;
box on;

subplot(2,2,2);
hold on;
plot(t1/3600, y1(:,2)*1e6);
axis([0,1,0,10]);
xlabel('Time (h)');
ylabel('[H_i] (\muM)');
axis square;
box on;

subplot(2,2,3);
hold on;
plot(t1/3600, y1(:,3)*1e6);
axis([0,1,0,7]);
xlabel('Time (h)');
ylabel('[KatG] (\muM)');
axis square;
box on;

subplot(2,2,4);
hold on;
plot(t1/3600, y1(:,4)*1e6);
axis([0,1,0,5]);
xlabel('Time (h)');
ylabel('[AhpCF] (\muM)');
axis square;
box on;
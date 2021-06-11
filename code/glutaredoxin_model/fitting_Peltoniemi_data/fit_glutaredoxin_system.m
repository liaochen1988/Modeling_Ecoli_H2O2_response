clear all;
clc;

%% Data from deglutationylation experiments
% Peltoniemi et al. (2006) The Journal of Biological Chemistry

% Experiment 1: GSH = 1 mM, Peptide = 5 uM, GSSG = 0 uM
Grx1 = [0,5,10,15,20,25]; %nM
initR1 = [0.725,9.425,17.578,26.845,35.980,45.121]; % nM/s

% Experiment 2: GSH = 1 mM, Grx1 = 20 nM, GSSG = 0 uM
Peptide2 = [0,1,2,3,4,5,6,8,10,15,20]; % uM
initR2 = [0.184,10.538,17.969,22.339,28.184,37.620,41.312,47.566,52.669,58.471,62.921];

% Experiment 3: Peptide = 5 uM, Grx1 = 20 nM, GSSG = 0 uM
Gsh3 = [0,0.05,0.1,0.25,0.5,1.0,1.5,2.0,3.0,4.0]; % mM
initR3 = [0.031,0.784,2.443,12.266,23.233,35.642,41.288,45.819,44.301,46.630]; % nM/s

% Experiment 4: Peptide = 5 uM, Grx1 = 20 nM, Gsh = 1 mM
Gssg4 = [0,5,10,15,25]; % uM
initR4 = [35.595,24.736,18.530,12.998,7.808]; % nM/s

%% Fit glutaredoxin model to experimental data
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [-1,0,0]; % k1+, k2+, k2-
[sol,~,~,exitflag] = lsqnonlin(@mse_grx_gsh,p0,[],[],options,Grx1,initR1,Peptide2,initR2,Gsh3,initR3,Gssg4,initR4);
assert(exitflag>0);

k1p = 10.^sol(1);
k2p = 10.^sol(2);
k2n = 10.^sol(3);

%% Plot fitting results

figure();

% Experiment 1: GSH = 1 mM, Peptide = 5 uM, GSSG = 0 uM
subplot(2,2,1);
hold on;

Grx = [0:0.1:25]*1e-3; % uM
Gsh = 1*1e3; % uM
Peptide = 5; % uM
Gssg = 0; % uM
ypred1 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

plot(Grx1,initR1,'ko');
plot(Grx*1e3,ypred1,'k-');
axis square;
box on;
axis([0,25,0,50]);
xlabel('Grx (nM)');
ylabel('Init Rate (nM/s)');

% Experiment 2: GSH = 1 mM, Grx1 = 20 nM, GSSG = 0 uM
subplot(2,2,2);
hold on;

Grx = 20*1e-3; % uM
Gsh = 1*1e3; % uM
Peptide = [0:0.1:20]; % uM
Gssg = 0; % uM
ypred2 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

plot(Peptide2,initR2,'ko');
plot(Peptide,ypred2,'k-');
axis square;
box on;
axis([0,20,0,80]);
xlabel('Peptide (\muM)');
ylabel('Init Rate (nM/s)');

% Experiment 3: Peptide = 5 uM, Grx1 = 20 nM, GSSG = 0 uM
subplot(2,2,3);
hold on;

Grx = 20*1e-3; % uM
Gsh = [0:0.1:4]*1e3; % uM
Peptide = 5; % uM
Gssg = 0; % uM
ypred3 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

plot(Gsh3,initR3,'ko');
plot(Gsh*1e-3,ypred3,'k-');
axis square;
box on;
axis([0,4,0,50]);
xlabel('GSH (mM)');
ylabel('Init Rate (nM/s)');

% Experiment 4: Peptide = 5 uM, Grx1 = 20 nM, GSH = 1 mM
subplot(2,2,4);
hold on;

Grx = 20*1e-3; % uM
Gsh = 1*1e3; % uM
Peptide = 5; % uM
Gssg = [0:0.1:25]; % uM
ypred4 = calc_init_rate(Grx,Gsh,Peptide,Gssg,k1p,k2p,k2n)*1e3; % nM/s

plot(Gssg4,initR4,'ko');
plot(Gssg,ypred4,'k-');
axis square;
box on;
axis([0,25,0,40]);
xlabel('GSSG (\muM)');
ylabel('Init Rate (nM/s)');

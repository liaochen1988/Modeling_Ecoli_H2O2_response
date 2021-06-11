clear all;
clc;

kon = 1.38e-5/60; % 1/uM/s
koff = 9.89e-6/60; % 1/uM/s

k1p = 0.51; % 1/uM/s
k1n = 0; % 1/uM^2/s
k2p = 5.15e-6; % 1/uM^2/s
k2n = 0.74; % 1/uM/s

MetE_SH = 50;  % uM
GSH = 5.0*1e3; % uM
GSSG = GSH/2; % uM
Grx = 0; % uM

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,4),'NonNegative',1:4);
tspan = [0,1e4]; % s

x0 = [MetE_SH,GSSG,0,GSH];
[t,y] = ode15s(@metE_activity_glutaredoxin,tspan,x0,options,kon,koff,k1p,k1n,k2p,k2n,Grx);

figure();

for i=1:4
    subplot(2,2,i);
    plot(t,y(:,i));
    axis square;
    box on;
end
clear all;
clc;

kon = 1.38e-5/60; % 1/uM/s
koff = 9.89e-6/60; % 1/uM/s

k1p = 0.51; % 1/uM/s
k1n = 0; % 1/uM^2/s
k2p = 5.15e-6; % 1/uM^2/s
k2n = 0.74; % 1/uM/s

MetE_SH = 5000; % uM
GSH = [1.1,2.3,4.6,9.1]*1e3; % uM
Grx = 20; % uM

GSH_GSSG_Ratio = 10.^[-2:0.1:2];
metE_activity = zeros(length(GSH),length(GSH_GSSG_Ratio));

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,4),'NonNegative',1:4);
tspan = [0,1e20]; % s

for i=1:length(GSH)
    for j=1:length(GSH_GSSG_Ratio)
        GSSG = GSH(i)/GSH_GSSG_Ratio(j);
        x0 = [MetE_SH,GSSG,0,GSH(i)];
        [t,y] = ode15s(@metE_activity_glutaredoxin,tspan,x0,options,kon,koff,k1p,k1n,k2p,k2n,Grx);
        metE_activity(i,j) = y(end,1)/MetE_SH;
    end
end

figure();
hold on;
cc = hsv(length(GSH));

for i=1:length(GSH)
    plot(GSH_GSSG_Ratio,metE_activity(i,:),'k-','Color',cc(i,:));
    plot(GSH_GSSG_Ratio,GSH_GSSG_Ratio./(1.4+GSH_GSSG_Ratio),'k--','Color',cc(i,:));
end

axis square;
box on;
set(gca,'XScale','log');
axis([1e-2,1e2,0,1]);
xlabel('[GSH]/[GSSG]');
ylabel('Relative MetE activity');
legend('1.1 mM GSH','2.3 mM GSH','4.6 mM GSH','9.1 mM GSH');
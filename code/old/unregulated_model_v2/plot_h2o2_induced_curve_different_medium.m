clear all;
clc;

% In vitro, glucose, glycerol, acetate
GSH  = [0.025, 0.0166,  0.0176,  0.00797]; % M
GSSG = [1e-4,  0.00237, 0.00731, 0.00168]; % M
cc = jet(length(GSH));

k_oxyr_ox = 9.9863; % 1/s
n_oxyr = 1.3614; 
Km_oxyr = 41.3179; % uM
k_oxyr_ox1 = 149.1049; % 1/s/M^4
k_oxyr_red1 = 1.5218e+10; % 1/s/M^8

figure();
hold on;

xdata_h2o2 = 10.^[-6:0.1:2]; % uM
for k=1:length(GSH)
    kswitch = k_oxyr_ox*xdata_h2o2.^n_oxyr./(Km_oxyr^n_oxyr+xdata_h2o2.^n_oxyr);
    ydata_foxyrox = (kswitch+k_oxyr_ox1*GSSG(k)^4)./(k_oxyr_red1*GSH(k)^8+kswitch+k_oxyr_ox1*GSSG(k)^4);
    plot(xdata_h2o2,ydata_foxyrox,'k-','Color',cc(k,:),'LineWidth',1);
end

set(gca,'XScale','log');
axis square;
box on;
xlabel('[H_2O_2] (\muM)');
ylabel('Oxidized OxyR fraction');
axis([1e-6,1e2,-0.1,1.1]);
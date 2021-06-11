clear all;
clc;

% rdec is decomposition rate of HPI
h2o2 = [0.5,  1.5,   3,     15,    30,     50,    150,    10000]; % uM
rdec = [-0.77,-0.71, -0.60, -0.47, -0.29,  -0.21, -0.074, -0.023]*(-1); % 1/min

% figure();
% plot(h2o2,rdec,'ko-');
% set(gca,'XScale','log');
% axis square;
% box on;
% xlabel('H2O2 (\muM)');
% ylabel('decomposition rate');

hill_fit = @(b,x)  b(1)./(1+(x/b(2)).^b(3));
b0 = [1,30,1];                                  % Initial Parameter Estimates
B = lsqcurvefit(hill_fit, b0, h2o2, rdec);
h2o2_fr = 10.^[-1:0.1:5];   % Plot Finer Resolution

% best fit
% B = [0.7696, 18.6285, 0.9483]

figure();
plot(h2o2, rdec, 'bp');
hold on;
plot(h2o2_fr, hill_fit(B,h2o2_fr), '-r');
hold off;
xlabel('H2O2 (\muM)');
ylabel('H2O2 decomposition rate (1/min)');
legend('Data', 'Hill Equation Fit', 'Location','NE');
set(gca,'XScale','log');

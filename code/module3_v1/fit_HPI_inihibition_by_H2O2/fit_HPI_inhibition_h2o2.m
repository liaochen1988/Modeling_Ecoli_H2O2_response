clear all;
clc;

% Experiment: Change of H2O2 decomposition rate at various H2O2
% SeaverImlay2001JB
% concentration in vitro (AhpCF dissociate and not functioning)
h2o2_seaver2001 = [0.5,  1.5,   3,     15,    30,     50,    150,    10000]; % uM
h2o2_deg_rate_seaver2001 = [-0.77,-0.71, -0.60, -0.47, -0.29,  -0.21, -0.074, -0.023]*(-1); % 1/min

hill_fit = @(b,x)  b(1)./(1+x/b(2));
b0 = [1,30]; % Initial Parameter Estimates
B = lsqcurvefit(hill_fit, b0, h2o2_seaver2001, h2o2_deg_rate_seaver2001);
h2o2_fr = 10.^[-1:0.1:5];   % Plot Finer Resolution

Ki_katG_h2o2 = B(2); % uM

figure();
plot(h2o2_seaver2001, h2o2_deg_rate_seaver2001, 'bp');
hold on;
plot(h2o2_fr, hill_fit(B,h2o2_fr), '-r');
hold off;
xlabel('H2O2 (\muM)');
ylabel('H2O2 decomposition rate (1/min)');
legend('Data', 'Hill Equation Fit', 'Location','NE');
set(gca,'XScale','log');

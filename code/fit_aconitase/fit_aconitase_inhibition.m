clear all;
clc;

h2o2 = [0 1 5 10 25 50 100 500]; % uM
activity = [100 93.994 74.749 38.796 13.554 4.698 3.463 2.321]/100;

hill_fit = @(b,x)  1./(1+(x/b(1)).^b(2));
b0 = [10,1]; % Initial Parameter Estimates
B = lsqcurvefit(hill_fit, b0, h2o2, activity);
h2o2_fr = 10.^[-1:0.1:5];   % Plot Finer Resolution

% best fit
% B = [8.3510 1.8042]

figure();
plot(h2o2, activity, 'bp');
hold on;
plot(h2o2_fr, hill_fit(B,h2o2_fr), '-r');
hold off;
xlabel('H2O2 (\muM)');
ylabel('Activity (%)');
legend('Data', 'Hill Equation Fit', 'Location','NE');
set(gca,'XScale','log');

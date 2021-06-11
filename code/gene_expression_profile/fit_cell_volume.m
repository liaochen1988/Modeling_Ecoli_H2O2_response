clear all;
clc;

growth_rate = [0.800, 0.806, 1.101, 1.376, 1.553, 1.840, 2.424] * log(2);
cell_volume = [0.424, 0.464, 0.589, 0.797, 1.039, 1.625, 2.736]; % um^3

hill_fit = @(b,x)  b(1)*exp(b(2)*x);
b0 = [1,1];    % Initial Parameter Estimates
B = lsqcurvefit(hill_fit, b0, growth_rate, cell_volume);
growth_rate_fg = [0:0.01:2];   % Plot Finer Resolution

% best fit
% B = [0.1882 1.6028]

figure();
plot(growth_rate, cell_volume, 'bp');
hold on;
plot(growth_rate_fg, hill_fit(B, growth_rate_fg), '-r');
hold off;
xlabel('Growth rate (h^{-1})');
ylabel('Cell volume (\mum^3)');
legend('Data', 'Exponential Fit', 'Location','NE');

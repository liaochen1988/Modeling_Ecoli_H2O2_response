clear all;
clc;

aeration = [0,25,32.5,50,100]/100;
h2o2_fmr_ext = [0,6.04,10.67,13.80,20.01]; % nM/min
h2o2_fmr_int = h2o2_fmr_ext*(1.2e-3)/(0.02*3.23e-15*5e8)/1e3/60;

%plot(aeration,h2o2_fmr_int,'ko-');

hill_fit = @(b,x)  b(1)*x./(b(2)+x);
b0 = [20,0.5];                                  % Initial Parameter Estimates
B = lsqcurvefit(hill_fit, b0, aeration, h2o2_fmr_int);
art_fr = [0:0.01:1];   % Plot Finer Resolution

% best fit
% B = [27.0252, 1.1546]

figure();
plot(aeration, h2o2_fmr_int, 'bp');
hold on;
plot(art_fr, hill_fit(B,art_fr), '-r');
hold off;
xlabel('aeration');
ylabel('H2O2 formation rate (\muM/s)');
legend('Data', 'Hill Equation Fit', 'Location','NE');

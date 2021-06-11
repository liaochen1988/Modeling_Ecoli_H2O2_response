clear all;
clc;

%% Experimental data
t_10 = [0, 10, 20, 30, 40, 50]; % min
H2O2_10 = [9.209, 3.914, 1.583, 0.579, 0.067, 0.028; ... % Wild type
           9.883, 7.302, 5.300, 3.443, 1.457, 0.512; ... % KatG
           9.838, 5.400, 2.263, 0.821, 0.154, 0.029; ... % KatE
           9.235, 7.087, 5.008, 3.074, 1.054, 0.232  ... % KatG and KatE
          ];
      
t_25 = [0, 10, 20, 30, 40, 50, 60]; % min
H2O2_25 = [25.803, 14.496, 8.445,  4.863,  2.394,  0.737,  0.130; ... % Wild type
           26.005, 23.503, 23.127, 22.759, 21.151, 19.942, 19.783; ... % KatG
           26.276, 15.196, 9.308, 5.350, 2.985, 1.088, 0.163; ... % KatE
           27.647, 23.142, 21.594, 23.385, 22.464, 22.331, 21.858 ... % KatG and KatE
           ];
       
t_100 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180];
H2O2_100 = [92.8512, 58.426, 38.009, 25.020, 14.808, 8.871, 5.303, 2.940, 0, 0, 0, 0, 0; ... % wild type
            92.109, 92.098, 86.274, 84.105, 80.022, 74.145, 70.291, 70.092, 65.375, 62.692, 59.144, 57.831, 56.645; ... % KatG
            92.932, 56.630, 35.869, 23.575, 14.923, 8.121, 4.714, 2.475, 0, 0, 0, 0, 0; ... % KatE
            92.298, 87.172, 86.132, 86.066, 84.835, 82.478, 81.553, 79.417, 79.578, 79.953, 78.124, 77.698, 74.130 ... % KatG and KatE
           ];

t_400 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240];
H2O2_400 = [377.918, 243.191, 171.930, 125.809, 97.431,  72.507,  61.637,  50.048,  36.496,  28.705, ...
            22.921, 17.524, 13.473, 8.308, 7.209, 5.190, 3.003; ... % wild type
            388.538, 289.409, 222.658, 194.062, 173.072, 151.092, 143.410, 125.679, 118.024, 108.138, ...
            97.458, 89.352, 85.853, 78.851, 74.656, 67.067, 63.623 ... % katE
           ];

%% Initial values of parameters to fit
alpha_max_katG  = 19.8/1e6/3600;        % M/s
alpha_max_ahpCF = 12.54/1e6/3600*1.5;   % M/s
alpha_max_grx   = 19.39/1e6/3600;       % M/s
vi_katG         = 2e-3;                 % 1/s
Ki_katG         = 4.2e-3;               % M
vi_ahpCF        = 2e-3;                 % 1/s
Ki_ahpCF        = 1e-5;                 % M

% %% Optimal values
% alpha_max_katG  = 5.5000e-09;           % M/s
% alpha_max_ahpCF = 5.2250e-09;           % M/s
% alpha_max_grx   = 5.3861e-09;           % M/s
% vi_katG         = 0.0117;               % 1/s
% Ki_katG         = 0.0041;               % M
% vi_ahpCF        = 0.0017;               % 1/s
% Ki_ahpCF        = 5.8779e-06;           % M

%% Fitting
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);

while i<=10
    [res,~,~,exitflag] = lsqnonlin(@mse_h2o2_challenge_v1, [vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF], [0,0,0,0], [Inf,Inf,Inf,Inf], options, alpha_max_katG, alpha_max_ahpCF, alpha_max_grx);
    assert(exitflag>0);
    
    vi_katG = res(1);
    Ki_katG = res(2);
    vi_ahpCF = res(3);
    Ki_ahpCF = res(4);
    
    [res,~,~,exitflag] = lsqnonlin(@mse_h2o2_challenge_v2, [alpha_max_katG, alpha_max_ahpCF, alpha_max_grx], [0,0,0], [Inf,Inf,Inf], options, vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF);
    assert(exitflag>0);
    
    alpha_max_katG = res(1);
    alpha_max_ahpCF = res(2);
    alpha_max_grx = res(3);
    
end

%% plot

[~, tspan, H2O2_output] = mse_h2o2_challenge_v1([vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF], alpha_max_katG, alpha_max_ahpCF, alpha_max_grx);

for i=1:size(H2O2_output,1)
    i
    subplot(2,2,i);
    hold on;
    plot(tspan / 3600, squeeze(H2O2_output(i,1,:)) * 1e6, '-r');
    plot(tspan / 3600, squeeze(H2O2_output(i,2,:)) * 1e6, '-b');
    
    axis square;
    box on;
    xlabel('Time (h)');
    ylabel('H_2O_2 (\muM)');
end

subplot(2,2,1);
plot(t_10 / 60, H2O2_10(1, :), 'ko', 'MarkerEdgeColor', 'r');
plot(t_10 / 60, H2O2_10(2, :), 'ko', 'MarkerEdgeColor', 'b');
axis([0, 1, 0, 12]);
title('10 \muM H_2O_2');

subplot(2,2,2);
plot(t_25 / 60, H2O2_25(1, :), 'ko', 'MarkerEdgeColor', 'r');
plot(t_25 / 60, H2O2_25(2, :), 'ko', 'MarkerEdgeColor', 'b');
axis([0, 1, 0, 35]);
title('25 \muM H_2O_2');

subplot(2,2,3);
plot(t_100 / 60, H2O2_100(1, :), 'ko', 'MarkerEdgeColor', 'r');
plot(t_100 / 60, H2O2_100(2, :), 'ko', 'MarkerEdgeColor', 'b');
axis([0, 3, 0, 120]);
title('100 \muM H_2O_2');

subplot(2,2,4);
plot(t_400 / 60, H2O2_400(1, :), 'ko', 'MarkerEdgeColor', 'r');
axis([0, 4, 0, 500]);
title('400 \muM H_2O_2');
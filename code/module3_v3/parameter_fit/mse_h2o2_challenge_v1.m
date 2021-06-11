function [dy, tspan, H2O2_output] = mse_h2o2_challenge_v1(p, alpha_max_katG, alpha_max_ahpCF, alpha_max_grx)

% Parameters
vi_katG = p(1);
Ki_katG = p(2);
vi_ahpCF = p(3);
Ki_ahpCF = p(4);

% Experimental data
t_10 = [0, 10, 20, 30, 40, 50]; % min
H2O2_10 = [9.209, 3.914, 1.583, 0.579, 0.067, 0.028; ... % Wild type
           9.883, 7.302, 5.300, 3.443, 1.457, 0.512; ... % KatG
           9.838, 5.400, 2.263, 0.821, 0.154, 0.029; ... % KatE
           9.235, 7.087, 5.008, 3.074, 1.054, 0.232  ... % KatG and KatE
          ]; % uM
      
t_25 = [0, 10, 20, 30, 40, 50, 60]; % min
H2O2_25 = [25.803, 14.496, 8.445,  4.863,  2.394,  0.737,  0.130; ... % Wild type
           26.005, 23.503, 23.127, 22.759, 21.151, 19.942, 19.783; ... % KatG
           26.276, 15.196, 9.308, 5.350, 2.985, 1.088, 0.163; ... % KatE
           27.647, 23.142, 21.594, 23.385, 22.464, 22.331, 21.858 ... % KatG and KatE
           ]; % uM
       
t_100 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]; % min
H2O2_100 = [92.8512, 58.426, 38.009, 25.020, 14.808, 8.871, 5.303, 2.940, 0, 0, 0, 0, 0; ... % wild type
            92.109, 92.098, 86.274, 84.105, 80.022, 74.145, 70.291, 70.092, 65.375, 62.692, 59.144, 57.831, 56.645; ... % KatG
            92.932, 56.630, 35.869, 23.575, 14.923, 8.121, 4.714, 2.475, 0, 0, 0, 0, 0; ... % KatE
            92.298, 87.172, 86.132, 86.066, 84.835, 82.478, 81.553, 79.417, 79.578, 79.953, 78.124, 77.698, 74.130 ... % KatG and KatE
           ]; % uM

t_400 = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180, 195, 210, 225, 240]; % min
H2O2_400 = [377.918, 243.191, 171.930, 125.809, 97.431,  72.507,  61.637,  50.048,  36.496,  28.705, ...
            22.921, 17.524, 13.473, 8.308, 7.209, 5.190, 3.003; ... % wild type
            388.538, 289.409, 222.658, 194.062, 173.072, 151.092, 143.410, 125.679, 118.024, 108.138, ...
            97.458, 89.352, 85.853, 78.851, 74.656, 67.067, 63.623 ... % katE
           ]; % uM

H2O2_IV = [9.209, 9.883, 9.838, 9.235;...
           25.803, 26.005, 26.276, 27.647;...
           92.8512, 92.109, 92.932, 92.298;...
           377.918, 400, 388.538, 400] * 1e-6; % M

%% simulation

H2O2_input = [10, 25, 100, 400] * 1e-6; % M
tspan = [0:0.1:4*3600];
H2O2_output = zeros(length(H2O2_input), 2, length(tspan));
N_cell_per_OD = 1e9; % cells/mL
Ve = 0.02; % L
InitN = 0.01 * N_cell_per_OD * Ve * 1e3;
options=odeset('RelTol',1e-10,'AbsTol',1e-10, 'NonNegative',[1:7]);

% Initial conditions
x0 = [0,0,0,0,0,0,InitN];
[~,x_wt] = ode15s(@h2o2_model, [0,1e10], x0, options, 0, false, alpha_max_katG, alpha_max_ahpCF, alpha_max_grx, vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF);
IC_wt = x_wt(end,:);

x0 = [0,0,0,0,0,0,InitN];
[~,x_katG] = ode15s(@h2o2_model, [0,1e10], x0, options, 0, true, alpha_max_katG, alpha_max_ahpCF, alpha_max_grx, vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF);
IC_katG = x_katG(end,:);

for i=1:length(H2O2_input)
    % wild type
    x0 = IC_wt;
    x0(1) = H2O2_IV(i,1);
    [~,x_wt] = ode15s(@h2o2_model, tspan, x0, options, 1, false, alpha_max_katG, alpha_max_ahpCF, alpha_max_grx, vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF);
    H2O2_output(i,1,:) = x_wt(:,1);
    
    % katG mutant
    x0 = IC_katG;
    x0(1) = H2O2_IV(i,2);
    [~,x_katG] = ode15s(@h2o2_model, tspan, x0, options, 1, true, alpha_max_katG, alpha_max_ahpCF, alpha_max_grx, vi_katG, Ki_katG, vi_ahpCF, Ki_ahpCF);
    H2O2_output(i,2,:) = x_katG(:,1);
end

dy = [];
% 10 uM H2O2
dy = [dy, (pchip(tspan, squeeze(H2O2_output(1,1,:)) * 1e6, t_10 * 60) - H2O2_10(1,:)) / max(H2O2_10(1,:))]; % WT
dy = [dy, (pchip(tspan, squeeze(H2O2_output(1,2,:)) * 1e6, t_10 * 60) - H2O2_10(2,:)) / max(H2O2_10(2,:))]; % KatG
% 25 uM H2O2
dy = [dy, (pchip(tspan, squeeze(H2O2_output(2,1,:)) * 1e6, t_25 * 60) - H2O2_25(1,:)) / max(H2O2_25(1,:))]; % WT
dy = [dy, (pchip(tspan, squeeze(H2O2_output(2,2,:)) * 1e6, t_25 * 60) - H2O2_25(2,:)) / max(H2O2_25(2,:))]; % KatG
% 100 uM H2O2
dy = [dy, (pchip(tspan, squeeze(H2O2_output(3,1,:)) * 1e6, t_100 * 60) - H2O2_100(1,:)) / max(H2O2_100(1,:))]; % WT
dy = [dy, (pchip(tspan, squeeze(H2O2_output(3,2,:)) * 1e6, t_100 * 60) - H2O2_100(2,:)) / max(H2O2_100(2,:))]; % KatG
% 400 uM H2O2
dy = [dy, (pchip(tspan, squeeze(H2O2_output(4,1,:)) * 1e6, t_400 * 60) - H2O2_400(1,:)) / max(H2O2_400(1,:))]; % WT

dy = dy';

end


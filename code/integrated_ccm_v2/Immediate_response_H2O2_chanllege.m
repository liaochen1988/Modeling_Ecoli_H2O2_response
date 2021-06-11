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

H2O2_IV = [9.209, 9.883, 9.838, 9.235;...
           25.803, 26.005, 26.276, 27.647;...
           92.8512, 92.109, 92.932, 92.298;...
           377.918, 400, 388.538, 400] * 1e-6;
       
%% Simulatioon options

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose (g/L)
glucose = 10 * 1e-3 * 180.156; % 10 mM

%   scaled dissolved oxygen level
a = 40; % aerobic condition

Num_of_State_Variable = 59;

V_env = 0.02; % L

%% Initial condition
% tolerance cannot be too tight
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);
IC = Initial_Concentration();
IC(2) = glucose;

% wild type
katG_mutant = false;
katE_mutant = false;
[~,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
IC_wt = x(end,:);

% katG mutant
katG_mutant = true;
katE_mutant = false;
[~,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC_wt,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
IC_katG = x(end,:);

% katE mutant
katG_mutant = false;
katE_mutant = true;
[~,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC_wt,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
IC_katE = x(end,:);

% katG and katE mutant
katG_mutant = true;
katE_mutant = true;
[~,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC_wt,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
IC_both = x(end,:);

%% simulation

H2O2_input = [10, 25, 100, 400] * 1e-6; % M
tspan = [0:1:4*3600];
H2O2_output = zeros(length(H2O2_input), 4, length(tspan));
N_cell_per_OD = 1e9;

options=odeset('RelTol',1e-10,'AbsTol',1e-10, 'NonNegative',[1:Num_of_State_Variable]);
for i=1:length(H2O2_input)
    i
    
    % wild type
    katG_mutant = false;
    katE_mutant = false;
    IC_wt(1) = 0.01 * N_cell_per_OD * V_env * 1e3;
    IC_wt(47) = H2O2_IV(i,1);
    [t_wt,x_wt] = ode15s(@Kinetic_model,tspan,IC_wt,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
    res = Kinetic_model(t_wt(end),x_wt(end,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false);
    H2O2_output(i,1,:) = x_wt(:, 47);
    
    % katG mutant
    katG_mutant = true;
    katE_mutant = false;
    IC_katG(1) = 0.01 * N_cell_per_OD * V_env * 1e3;
    IC_katG(47) = H2O2_IV(i,2);
    [t_katG,x_katG] = ode15s(@Kinetic_model,tspan,IC_katG,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
    res = Kinetic_model(t_katG(end),x_katG(end,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false);
    res(1) * 3600
    H2O2_output(i,2,:) = x_katG(:, 47);
    
    % katE mutant
    katG_mutant = false;
    katE_mutant = true;
    IC_katE(1) = 0.01 * N_cell_per_OD * V_env * 1e3;
    IC_katE(47) = H2O2_IV(i,3);
    [t_katE,x_katE] = ode15s(@Kinetic_model,tspan,IC_katE,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
    res = Kinetic_model(t_katE(end),x_katE(end,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false);
    res(1) * 3600
    H2O2_output(i,3,:) = x_katE(:, 47);
    
    % katE mutant
    katG_mutant = true;
    katE_mutant = true;
    IC_both(1) = 0.01 * N_cell_per_OD * V_env * 1e3;
    IC_both(47) = H2O2_IV(i,4);
    [t_both,x_both] = ode15s(@Kinetic_model,tspan,IC_both,options,arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,true);
    res = Kinetic_model(t_both(end),x_both(end,:),arcA_mutant,fnr_mutant,katE_mutant,katG_mutant,a,V_env,false);
    res(1) * 3600
    H2O2_output(i,4,:) = x_both(:, 47);
end
disp('Done');

%% plot 
figure();

cc = [0, 0, 0; ...
      30, 144, 255; ...
      34, 139, 34; ...
      138, 43, 226]/255;

% plot simulation
for i=1:length(H2O2_input)
    
    subplot(2,2,i);
    hold on;
    
    for j=1:4
        plot(tspan / 3600, squeeze(H2O2_output(i,j,:)) * 1e6, '-k', 'Color', cc(j,:));
    end
    
    axis square;
    box on;
    xlabel('Time (h)');
    ylabel('H_2O_2 (\muM)');
end

% plot experimental data
for i=1:4
    subplot(2,2,1);
    plot(t_10 / 60, H2O2_10(i, :), 'ko', 'MarkerEdgeColor', cc(i,:));
    subplot(2,2,2);
    plot(t_25 / 60, H2O2_25(i, :), 'ko', 'MarkerEdgeColor', cc(i,:));
    subplot(2,2,3);
    plot(t_100 / 60, H2O2_100(i, :), 'ko', 'MarkerEdgeColor', cc(i,:));
end
    
subplot(2,2,4);
plot(t_400 / 60, H2O2_400(1, :), 'ko', 'MarkerEdgeColor', cc(1,:));
plot(t_400 / 60, H2O2_400(2, :), 'ko', 'MarkerEdgeColor', cc(3,:));

subplot(2,2,1);
axis([0, 1, 0, 12]);
title('10 \muM H_2O_2');

subplot(2,2,2);
axis([0, 1, 0, 35]);
title('25 \muM H_2O_2');

subplot(2,2,3);
axis([0, 3, 0, 120]);
title('100 \muM H_2O_2');

subplot(2,2,4);
axis([0, 4, 0, 500]);
title('400 \muM H_2O_2');

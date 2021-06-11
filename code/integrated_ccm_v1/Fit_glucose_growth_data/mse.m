function [dx, Glucose_obs, Growth_rate_obs, GLC_ex_input, Growth_rate_sim] = mse(p)

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose 10 g/L
Glucose_obs =      [0.008, 0.011, 0.015, 0.038, 0.056, 0.070, 0.057, 0.049, 0.085, 0.123, ...
    0.170, 0.155, 0.170, 0.211, 0.285, 0.307, 0.372, 0.397, 0.496, 0.629, ...
    0.644, 0.336, 0.476, 0.550, 0.749, 0.838] * 1e-3 * 180.156; % mM to g/L
Growth_rate_obs =   [0.026, 0.078, 0.135, 0.283, 0.373, 0.398, 0.445, 0.535, 0.780, 0.850, ...
    0.863, 0.794, 0.814, 0.851, 0.817, 0.814, 0.818, 0.808, 0.820, 0.832, ...
    0.854, 0.959, 0.914, 0.897, 0.921, 0.871]; % 1/h

GLC_ex_input = 10.^[-4:0.2:0];
Growth_rate_sim = zeros(length(GLC_ex_input),1);

%   scaled dissolved oxygen level
a = 40;

integrated_model = true;
Num_of_State_Variable = 40;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

IC = Initial_Concentration();
for i=1:length(GLC_ex_input)
    IC(2) = GLC_ex_input(i);
    called_by_ode15s = true;
    [~,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,p);
    
    called_by_ode15s = false;
    Growth_rate_sim(i) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,p) * 3600;
    
    IC = x(end,:);
end

dx = pchip(GLC_ex_input, Growth_rate_sim, Glucose_obs) - Growth_rate_obs;

end


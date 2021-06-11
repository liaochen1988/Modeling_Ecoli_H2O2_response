clear all;
clc;

%% Fit k_ahp and Km_ahp
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-12,'TolFun',1e-12,'MaxFunEval',Inf,'Algorithm','levenberg-marquardt');
p0 = [170]; 
[c_met_e,~,~,exitflag] = lsqnonlin(@mse,p0,[0],[Inf],options);
assert(exitflag>0);
[~, Glucose_obs, Growth_rate_obs, GLC_ex_input, Growth_rate_sim] = mse(c_met_e);

%%  plot
figure();
hold on;

plot(Glucose_obs, Growth_rate_obs, 'ro');
plot(GLC_ex_input, Growth_rate_sim, 'k-');

axis square;
box on;
xlabel('Glucose (g/L)');
ylabel('Growth rate (1/h)');

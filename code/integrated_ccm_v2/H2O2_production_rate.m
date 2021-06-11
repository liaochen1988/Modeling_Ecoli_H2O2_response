clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose (g/L)
GLC_ex_input = 10;

%   scaled dissolved oxygen level
a = [1:4, 5:5:100]; % aerobic condition

growthRate    = zeros(length(a),1);
H2O2_prod_rate = zeros(length(a),1);

%% simulation
Num_of_State_Variable = 59;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

IC = Initial_Concentration();
for i=1:length(a)
    [i]
    IC(2) = GLC_ex_input;
    called_by_ode15s = true;
    [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a(i),called_by_ode15s);
    
    called_by_ode15s = false;
    res = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a(i),called_by_ode15s);
    growthRate(i) = res(1);
    H2O2_prod_rate(i) = res(2);
    
    IC = x(end,:);
end

%% plot growth rate
figure();

subplot(1,2,1);
hold on;
plot(a, growthRate * 3600);
axis square;
box on;
ylabel('Growth rate (1/h)');
xlabel('Aerobiosis');

subplot(1,2,2);
hold on;
plot(a, H2O2_prod_rate * 1e6);
axis square;
ylabel('H_2O_2 prod. rate (\muM/s)');
xlabel('Aerobiosis');
box on;



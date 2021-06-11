clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose 10 g/L
GLC_ex_input = [10.^[-5:0.5:0]];

%   scaled dissolved oxygen level
a = 1;

Num_of_State_Variable = 36;
Num_of_Flux = 11;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

growth_rate = zeros(length(GLC_ex_input),1);
IC = Initial_Concentration();
for i=1:length(GLC_ex_input)
    i
    IC(2) = GLC_ex_input(i);
    called_by_ode15s = true;
    [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,called_by_ode15s);
    
    called_by_ode15s = false;
    growth_rate(i) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a,called_by_ode15s);
    
    IC = x(end,:);
end

%%  plot
figure();
hold on;
plot(GLC_ex_input, growth_rate,'ko-');
axis square;
box on;

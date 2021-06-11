clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose 10 g/L
GLC_ex_input = 10.^[-4:0.5:0];

%   scaled dissolved oxygen level
a = [0.2:0.2:0.8,1.0:1.0:15.0,20.0:10.0:50.0];

integrated_model = true;
Num_of_State_Variable = 40;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

oxygen_consumption_rate = zeros(length(GLC_ex_input),length(a));

c_met_e = 170;
IC = Initial_Concentration();
for i=1:length(GLC_ex_input)
    IC(2) = GLC_ex_input(i);
    gr_temp = zeros(1,length(a));
    for j=1:length(a)   
        [i,j]
        
        called_by_ode15s = true;
        [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a(j),called_by_ode15s,integrated_model,c_met_e);
        
        called_by_ode15s = false;
        oxygen_consumption_rate(i,j) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a(j),called_by_ode15s,integrated_model,c_met_e);
        
        IC = x(end,:);
    end
end

%%  plot
figure();
hold on;

cc = jet(length(GLC_ex_input));
for i=1:length(GLC_ex_input)
    plot(a, oxygen_consumption_rate(i,:), 'k-o', 'Color', cc(i,:));
end

axis square;
box on;
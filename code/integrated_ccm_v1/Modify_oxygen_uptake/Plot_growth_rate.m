clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose 10 g/L
GLC_ex_input = [1e-3]; %10.^[-5:0.5:0];

%   scaled dissolved oxygen level
a = [1,100];%10.^[0:0.5:2];

growth_rate = zeros(length(a),length(GLC_ex_input),2);

c_met_e = 170;

%% original oxygen uptake rate
new_qOCR = false;
Num_of_State_Variable = 40;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

for k=1:length(a)
    IC = Initial_Concentration();
    IC = IC(1:40);
    for i=1:length(GLC_ex_input)
        [k,i]
        IC(2) = GLC_ex_input(i);
        called_by_ode15s = true;
        [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a(k),called_by_ode15s,new_qOCR,c_met_e);
        
        called_by_ode15s = false;
        growth_rate(k,i,1) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a(k),called_by_ode15s,new_qOCR,c_met_e);
        
        IC = x(end,:);
    end
end

%% modified oxygen uptake rate
new_qOCR = true;
Num_of_State_Variable = 41;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

for k=1:length(a)
    IC = Initial_Concentration();
    for i=1:length(GLC_ex_input)
        [k i]
        IC(2) = GLC_ex_input(i);
        called_by_ode15s = true;
        [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a(k),called_by_ode15s,new_qOCR,c_met_e);
        
        called_by_ode15s = false;
        growth_rate(k,i,2) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a(k),called_by_ode15s,new_qOCR,c_met_e);
        
        IC = x(end,:);
    end
end

%%  plot
figure();
hold on;

cc = jet(length(a));
for k=1:length(a)
    plot(GLC_ex_input, squeeze(growth_rate(k,:,1)), 'ko-', 'MarkerEdgeColor', cc(k,:));
    plot(GLC_ex_input, squeeze(growth_rate(k,:,2)), 'kx-', 'MarkerEdgeColor', cc(k,:));
end
axis square;
box on;
xlabel('Glucose (g/L)');
ylabel('Growth rate (1/h)');

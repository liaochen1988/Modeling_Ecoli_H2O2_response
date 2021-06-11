clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose 10 g/L
GLC_ex_input = [10.^[-5:0.5:0]];
growth_rate_ori = zeros(length(GLC_ex_input),1);

%   scaled dissolved oxygen level
a = 1;

%% original model
integrated_model = false;
Num_of_State_Variable = 36;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

IC = Initial_Concentration();
IC = IC(1:36);
for i=1:length(GLC_ex_input)
    i
    IC(2) = GLC_ex_input(i);
    called_by_ode15s = true;
    [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,0);
    
    called_by_ode15s = false;
    growth_rate_ori(i,1) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,0);
    
    IC = x(end,:);
end

%% integrated model
integrated_model = true;
Num_of_State_Variable = 39;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

c_met_e = 10.^[2:0.01:2.2];
growth_rate_int = zeros(length(GLC_ex_input), length(c_met_e));

for k=1:length(c_met_e)
    IC = Initial_Concentration();
    for i=1:length(GLC_ex_input)
        [k i]
        IC(2) = GLC_ex_input(i);
        called_by_ode15s = true;
        [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,c_met_e(k));
        
        called_by_ode15s = false;
        growth_rate_int(i,k) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,c_met_e(k));
        
        IC = x(end,:);
    end
end

%%  plot
figure();
hold on;
plot(GLC_ex_input, growth_rate_ori,'ko-');

cc = jet(length(c_met_e));
for i=1:length(c_met_e)
    plot(GLC_ex_input, growth_rate_int(:,i),'ko-','MarkerFaceColor',cc(i,:));
end
axis square;
box on;
xlabel('Glucose (g/L)');
ylabel('Growth rate (1/h)');

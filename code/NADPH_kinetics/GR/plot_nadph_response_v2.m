clear all;
clc;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose (g/L)
glucose = 2;

%   scaled dissolved oxygen level
a = 40; % aerobic condition

Num_of_State_Variable = 37;
tspan = [0:0.1:60];

%% Running to steady state
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);
IC = Initial_Concentration();
IC(2) = glucose;

[t,x] = ode15s(@Kinetic_model,[-10000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,true,0,0);
IC = x(end,:);

%% Run NADPH response
%GSSG = 0.00237; % M
GSSG = [5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2];     % M
% Gor = 2.8499 * 1e-6;         % M
Gor = [5e-7, 1e-6, 5e-6, 1e-5];
cls = jet(length(GSSG));

figure();

for i=1:length(Gor)
    for j=1:length(GSSG)
        [i,j]
        [t,x] = ode15s(@Kinetic_model,[0.0 60.0],IC,options,arcA_mutant,fnr_mutant,a,true,GSSG(j),Gor(i));
        
        subplot(2,4,i);
        hold on;
        plot(t,x(:,37) * 564,'k-','Color',cls(j,:));
        
        subplot(2,4,i+4);
        hold on;
        plot(t,x(:,37)/x(1,37),'k-','Color',cls(j,:));
    end
    
    subplot(2,4,i);
    hold on;
    axis([-1,60,0,150]);
    plot([-1,60],[123,123],'k--');
    xlabel('Time (s)');
    ylabel('NADPH (uM)');
    axis square;
    box on;
    title(sprintf('Gor: %2.2f uM', Gor(i)*1e6));
    
    subplot(2,4,i+4);
    hold on;
    
    axis([-1,60,0,1.5]);
    xlabel('Time (s)');
    ylabel('Relative NADPH');
    
    t_1mM = [0, 5, 15, 30, 60]; % sec
    NADPH_1mM = [1.000	1.163	0.645	0.697	0.568]; % Relative
    plot(t_1mM, NADPH_1mM, 'ko');
    axis square;
    box on;
    title(sprintf('Gor: %2.2f uM', Gor(i)*1e6));
end
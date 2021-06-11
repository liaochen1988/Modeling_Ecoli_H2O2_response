clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose (g/L)
glucose = 10;

%   scaled dissolved oxygen level
a = 40; % aerobic condition

Num_of_State_Variable = 59;


%% presimulation
% tolerance cannot be too tight
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);
IC = Initial_Concentration();
IC(2) = glucose;
[~,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,true);
IC = x(end,:);

%% simulation
Hout = [0, 10.^[-9:0.5:0]]; % M
Hin_fw = zeros(size(Hout));
growthRate = zeros(size(Hout));
OxyR_ox = zeros(size(Hout));
OxyR_red = zeros(size(Hout));

options=odeset('RelTol',1e-10,'AbsTol',1e-10, 'NonNegative',[1:Num_of_State_Variable]);
for i=1:length(Hout)
    i
    IC(47) = Hout(i);
    [t,x] = ode15s(@Kinetic_model,[-10000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,true);
    res = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a,false);
    growthRate(i) = res(1);
    Hin_fw(i) = x(end,48);
    OxyR_ox(i) = x(end,49);
    OxyR_red(i) = x(end,50);
    
%     res(1)
%     res(2)
%     Hin_fw(i)
%     OxyR_ox(i)/(OxyR_ox(i)+OxyR_red(i))
%     
    IC = x(end,:);
end

%% plot growth rate
%figure();

subplot(2,2,1);
hold on;
plot(Hout * 1e6, growthRate * 3600);
set(gca,'XScale','log');
axis square;
box on;
ylabel('Growth rate (1/h)');
xlabel('Extracellular H_2O_2 (\muM)');

subplot(2,2,2);
hold on;
plot(Hout * 1e6, Hin_fw * 1e6);
set(gca,'XScale','log');
set(gca,'YScale','log');
axis square;
ylabel('Intracellular H_2O_2 (\muM)');
xlabel('Extracellular H_2O_2 (\muM)');
box on;

subplot(2,2,3);
hold on;
plot(Hout * 1e6, OxyR_ox * 1e6, 'r-');
plot(Hout * 1e6, OxyR_red * 1e6, 'b-');
plot(Hout * 1e6, (OxyR_ox + OxyR_red) * 1e6, 'k-');
set(gca,'XScale','log');
axis square;
ylabel('Oxidized OxyR (\muM)');
xlabel('Extracellular H_2O_2 (\muM)');
box on;

subplot(2,2,4);
hold on;
plot(Hout * 1e6, OxyR_ox ./ (OxyR_ox + OxyR_red), 'r-');
set(gca,'XScale','log');
axis square;
ylabel('Fraction of oxidized OxyR');
xlabel('Extracellular H_2O_2 (\muM)');
box on;
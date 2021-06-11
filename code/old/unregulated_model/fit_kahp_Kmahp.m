clear all;
clc;

%% best fit parameter values
% k_ahp = 1.1154e-18;
% Km_ahp = 5.3148e-10;

%% parameters
P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20;
k_cat = 2.7e-13; % L/s
k_ahp = 2.1e-18; % mol/s
Km_ahp = 1.2e-6; % M
n_ahp = 1;

%% oxyR oxidation data
h2o2_invitro = [0,0.025,0.05,0.075,0.1,0.2,0.5,1.0,2.0,5.0,10.0]/1e6; % M
oxyr_ox_invitro = [0,0.216,0.300,0.467,0.533,0.732,0.865,0.937,0.966,0.978,1.000];

h2o2_invivo = [0,0.5,1.0,2.0,5.0,10.0]/1e6; % M
oxyr_ox_invivo = [0,0,0,0,0.6,1.0];

%% fit Hill function
options = optimoptions('lsqnonlin','Display','off','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [1.5,0.1];
lb = [0,0];
ub = [Inf,Inf];
[hill_opt_invitro,~,~,exitflag] = lsqnonlin(@mseHill,p0,lb,ub,options,h2o2_invitro,oxyr_ox_invitro);
assert(exitflag>0);

p0 = [10.7,4.8e-6];
lb = [0,0];
ub = [Inf,Inf];
[hill_opt_invivo,~,~,exitflag] = lsqnonlin(@mseHill,p0,lb,ub,options,h2o2_invivo,oxyr_ox_invivo);
assert(exitflag>0);

%% Fit k_ahp and Km_ahp
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf,'Algorithm','levenberg-marquardt');
p0 = [log10(k_ahp),log10(Km_ahp)];
lb = [-Inf,-Inf];
ub = [Inf,Inf];
[p_opt,~,~,exitflag] = lsqnonlin(@mseOxyR,p0,lb,ub,options,k_met,PA,k_cat,n_ahp,hill_opt_invitro,h2o2_invivo,oxyr_ox_invivo);
assert(exitflag>0);
k_ahp = 10.^(p_opt(1));
Km_ahp = 10.^(p_opt(2));

%% Predict OxyR oxidation percentage
Hout = 10.^[-9:0.001:-4]; % M
Hin = zeros(size(Hout));
options = optimset('Display','off','TolX',1e-12); % show iterations
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin(i+1);
    end
    [Hin(i),fval,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,k_ahp,Km_ahp,n_ahp);
    while Hin(i)<0
        x0 = 10.^(-9+rand*5);
        [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,k_ahp,Km_ahp,n_ahp);
    end
    assert(exitflag>0);
end
OxyR = Hin.^hill_opt_invitro(1)./(Hin.^hill_opt_invitro(1)+hill_opt_invitro(2)^hill_opt_invitro(1));

%% plot Hout vs OxyR
figure();
hold on;
plot(h2o2_invitro*1e6, oxyr_ox_invitro, 'ko','MarkerFaceColor','b');
plot(Hin*1e6, OxyR,'b-');
plot(h2o2_invivo*1e6, oxyr_ox_invivo, 'ko','MarkerFaceColor','r');
plot(Hout*1e6, OxyR, 'r-');
legend('In vitro obs','In vitro sim','In vivo obs','In vivo sim');
title('Simulation');
axis([1e-3,1e2,-0.1,1.1]);
axis square;
box on;
ylabel('Oxidized OxyR (%)');
xlabel('Extracellular H_2O_2 (\muM)');
set(gca,'XScale','log');
set(gca,'XTick',10.^[-3:1:2]);
set(gca,'YTick',[0:0.2:1.0]);
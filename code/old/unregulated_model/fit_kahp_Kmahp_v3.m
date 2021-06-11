clear all;
clc;

%% best fit parameter values
% k_ahp = 1.1213e-18;
% Km_ahp = 6.2611e-10;

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

%% Fit k_ahp and Km_ahp
n_oxyr = 1.3614;
Km_oxyr = 41.3179/1e6; % uM
kox_oxyr= 9.9863; % 1/s
kred_oxyr = 0.0023; % 1/s, corresponding to half life 5 min
hill_opt_invitro = [n_oxyr,Km_oxyr,kox_oxyr,kred_oxyr];
    
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-12,'TolFun',1e-12,'MaxFunEval',Inf,'Algorithm','levenberg-marquardt');
p0 = [log10(k_ahp),log10(Km_ahp)];
lb = [-Inf,-Inf];
ub = [0,0];
[p_opt,~,~,exitflag] = lsqnonlin(@mseOxyR_v2,p0,lb,ub,options,k_met,PA,k_cat,n_ahp,hill_opt_invitro,h2o2_invivo,oxyr_ox_invivo);
assert(exitflag>0);
k_ahp = 10.^(p_opt(1));
Km_ahp = 10.^(p_opt(2));

%% change Km_ahp
Km_ahp_arr = 10.^[-9:1:-4];
Hout = 10.^[-9:0.001:-4]; % M
cc = jet(length(Km_ahp_arr));
options = optimset('Display','off','TolX',1e-12); % show iterations

for k=1:length(Km_ahp_arr)
    Hin = zeros(size(Hout));
    for i=length(Hout):-1:1
        if(i==length(Hout))
            x0 = Hout(i);
        else
            x0 = Hin(i+1);
        end
        [Hin(i),fval,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,k_ahp,Km_ahp_arr(k),n_ahp);
        while Hin(i)<0
            x0 = 10.^(-9+rand*5);
            [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,k_ahp,Km_ahp_arr(k),n_ahp);
        end
        assert(exitflag>0);
    end
    oxyr_ox_pred = Hin.^n_oxyr./(Hin.^n_oxyr*(1+kred_oxyr/kox_oxyr)+kred_oxyr/kox_oxyr*Km_oxyr^n_oxyr);
    
%     subplot(1,3,1);
%     hold on;
%     plot(Hout,Hin./Hout,'k-','Color',cc(k,:));
%     axis square;
%     box on;

    subplot(1,2,1);
    hold on;
    plot(Hout*1e6,oxyr_ox_pred,'k-','Color',cc(k,:));
    axis square;
    box on;
end

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
oxyr_ox_pred = Hin.^n_oxyr./(Hin.^n_oxyr*(1+kred_oxyr/kox_oxyr)+kred_oxyr/kox_oxyr*Km_oxyr^n_oxyr);

%% plot Hout vs OxyR
subplot(1,2,2);
hold on;
plot(h2o2_invitro*1e6, oxyr_ox_invitro, 'ko','MarkerFaceColor','b');
plot(Hin*1e6, oxyr_ox_pred,'b-');
plot(h2o2_invivo*1e6, oxyr_ox_invivo, 'ko','MarkerFaceColor','r');
plot(Hout*1e6, oxyr_ox_pred, 'r-');
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
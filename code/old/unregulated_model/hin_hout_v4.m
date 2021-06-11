clear all;
clc;

P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20; % mol/s, corresponding to 14 uM/s
%k_met = 0;
Km_ahp = 1.4138e-07; % M

Hout = 10.^[-12:0.1:0]; % M

% wild type
k_cat = 2.7e-13; % L/s
k_ahp = 2.1e-18; % mol/s
Delta = Km_ahp*PA+Km_ahp*k_cat-Hout*PA-k_met+k_ahp;
Hin_wt = (-Delta+sqrt(Delta.^2+4*Km_ahp*(PA+k_cat)*(k_met+Hout*PA)))/2/(PA+k_cat);

% catalase mutant
k_cat = 0;
k_ahp = 2.1e-18; % mol/s
Delta = Km_ahp*PA+Km_ahp*k_cat-Hout*PA-k_met+k_ahp;
Hin_cat_ = (-Delta+sqrt(Delta.^2+4*Km_ahp*(PA+k_cat)*(k_met+Hout*PA)))/2/(PA+k_cat);

% Ahp mutant
k_cat = 2.7e-13;
k_ahp = 0;
Delta = Km_ahp*PA+Km_ahp*k_cat-Hout*PA-k_met+k_ahp;
Hin_ahp_ = (-Delta+sqrt(Delta.^2+4*Km_ahp*(PA+k_cat)*(k_met+Hout*PA)))/2/(PA+k_cat);

% Catalase and Ahp double mutant
k_cat = 0;
k_ahp = 0;
Delta = Km_ahp*PA+Km_ahp*k_cat-Hout*PA-k_met+k_ahp;
Hin_cat_ahp_ = (-Delta+sqrt(Delta.^2+4*Km_ahp*(PA+k_cat)*(k_met+Hout*PA)))/2/(PA+k_cat);

% plot Hout vs Hin
%figure();
subplot(1,2,2);
hold on;

plot(Hout*1e6, Hin_wt*1e6, 'k-');
plot(Hout*1e6, Hin_cat_*1e6, 'r-');
plot(Hout*1e6, Hin_ahp_*1e6, 'b-');
plot(Hout*1e6, Hin_cat_ahp_*1e6, 'g-');
axis([1e-3,1e6,1e-2,1e6]);
axis square;
box on;
ylabel('Intracellular H_2O_2 (\muM)');
xlabel('Extracellular H_2O_2 (\muM)');
set(gca,'XScale','log');
set(gca,'YScale','log');
legend('WT','Cat-','Ahp-','Cat-Ahp-');
set(gca,'XTick',10.^[-3:3:6]);
set(gca,'YTick',10.^[-3:3:6]);

%% plot Hout vs Fold change
%figure();
subplot(1,2,1);
hold on;

plot(Hout*1e6, Hin_wt./Hout, 'k-');
plot(Hout*1e6, Hin_cat_./Hout, 'r-');
plot(Hout*1e6, Hin_ahp_./Hout, 'b-');
plot(Hout*1e6, Hin_cat_ahp_./Hout, 'g-');
axis([1e-3,1e6,0,1]);
axis square;
box on;
ylabel('Fold change');
xlabel('Extracellular H_2O_2 (\muM)');
set(gca,'XScale','log');
legend('WT','Cat-','Ahp-','Cat-Ahp-');
set(gca,'XTick',10.^[-3:3:6]);
set(gca,'YTick',[0:0.2:1]);
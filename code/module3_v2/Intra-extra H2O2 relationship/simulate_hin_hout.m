clear all;
clc;

P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20; % mol/s
Vc = 3.23*1e-15; % L
Ve = 0.2e-3;  % L

k_kat = 1.63e4; % 1/s
Km_kat = 3.9/1e3;% M
Ki_h2o2 = 19.55/1e6; %M
alpha_max_kat = 19.8/1e6; % M/h

k_ahp = 52.4; % 1/s
Km_ahp = 1.4987e-07; % M
alpha_max_ahp = 12.54/1e6; % M/h

alpha_max_grx = 19.39/1e6; % M/h

n_oxyr = 1.3614;
Km_oxyr = 41.3179e-6; % M
kon_oxyr_h2o2 = 9.99; % 1/s
kon_oxyr_grx = 1.49e7; % 1/s/M^5
koff_oxyr_grx = 1.52e15; % 1/s/M^9

GSSG = 1e-4; % M
GSH = 0.025; % M
Grx = 10*1e-6; % M
lambda = 0.99/3600; % 1/s
N = 0; %0.43*5e8; 

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,5),'NonNegative',1:5);

%% Forward simulation

Hout = 10.^[-9:0.1:-3,-3:0.01:0,0:0.1:3]; % M
Hin_fw = zeros(size(Hout));
for i=1:length(Hout)
    i
    if (i==1)
        x0 = [Hout(i),1,1,1,1];
    else
        x0 = y(end,:);
        x0(1) = Hout(i);
    end
    tspan = [0,1e10];
    [~,y] = ode15s(@h2o2_model,tspan,x0,options,N,GSH,GSSG,lambda,...
                   PA,k_met,Vc,Ve,...
                   k_kat,Km_kat,Ki_h2o2,alpha_max_kat,...
                   k_ahp,Km_ahp,alpha_max_ahp,...
                   alpha_max_grx,...
                   n_oxyr,Km_oxyr,kon_oxyr_h2o2,kon_oxyr_grx,koff_oxyr_grx);
    Hin_fw(i) = y(end,2);
end

kswitch = kon_oxyr_h2o2*Hin_fw.^n_oxyr./(Km_oxyr^n_oxyr+Hin_fw.^n_oxyr);
f_fw = (kswitch+kon_oxyr_grx*Grx*GSSG.^4)./(kswitch+koff_oxyr_grx*Grx*GSH.^8+kon_oxyr_grx*Grx*GSSG.^4);

%% Forward simulation

Hout = 10.^[-9:0.1:-3,-3:0.01:0,0:0.1:3]; % M
Hin_rv = zeros(size(Hout));
for i=length(Hout):-1:1
    i
    if (i==length(Hout))
        x0 = [Hout(i),1,1,1,1];
    else
        x0 = y(end,:);
        x0(1) = Hout(i);
    end
    
    tspan = [0,1e10];
    [~,y] = ode15s(@h2o2_model,tspan,x0,options,N,GSH,GSSG,lambda,...
                   PA,k_met,Vc,Ve,...
                   k_kat,Km_kat,Ki_h2o2,alpha_max_kat,...
                   k_ahp,Km_ahp,alpha_max_ahp,...
                   alpha_max_grx,...
                   n_oxyr,Km_oxyr,kon_oxyr_h2o2,kon_oxyr_grx,koff_oxyr_grx);
    Hin_rv(i) = y(end,2);
end

kswitch = kon_oxyr_h2o2*Hin_rv.^n_oxyr./(Km_oxyr^n_oxyr+Hin_rv.^n_oxyr);
f_rv = (kswitch+kon_oxyr_grx*Grx*GSSG.^4)./(kswitch+koff_oxyr_grx*Grx*GSH.^8+kon_oxyr_grx*Grx*GSSG.^4);

%% plot

figure();
hold on;

yyaxis left
plot(Hout*1e6,Hin_fw*1e6);
plot(Hout*1e6,Hin_rv*1e6);
xlabel('[H_o] (\muM)');
ylabel('[H_i] (\muM)');
set(gca,'XScale','log');
set(gca,'YScale','log');

yyaxis right
plot(Hout*1e6,f_fw);
plot(Hout*1e6,f_rv);
ylabel('Fraction of oxidized OxyR');

axis square;
box on;


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

%% Run forward and backward simulations

alpha_max_kat = [0:1:30]/1e6; % M/h
Hout = 10.^[-9:0.1:3]; % M
Hin_fw = zeros(length(alpha_max_kat),length(Hout));
Hin_rv = zeros(length(alpha_max_kat),length(Hout));

%% Forward simulation
for k=1:length(alpha_max_kat)
    k
    for i=1:length(Hout)
        if (i==1)
            x0 = [Hout(i),1,1,1,1];
        else
            x0 = y(end,:);
            x0(1) = Hout(i);
        end
        tspan = [0,1e10];
        [~,y] = ode15s(@h2o2_model,tspan,x0,options,N,GSH,GSSG,lambda,...
            PA,k_met,Vc,Ve,...
            k_kat,Km_kat,Ki_h2o2,alpha_max_kat(k),...
            k_ahp,Km_ahp,alpha_max_ahp,...
            alpha_max_grx,...
            n_oxyr,Km_oxyr,kon_oxyr_h2o2,kon_oxyr_grx,koff_oxyr_grx);
        Hin_fw(k,i) = y(end,2);
    end
end

%% Reverse simulation
for k=1:length(alpha_max_kat)
    k
    for i=length(Hout):-1:1
        if (i==length(Hout))
            x0 = [Hout(i),1,1,1,1];
        else
            x0 = y(end,:);
            x0(1) = Hout(i);
        end
        
        tspan = [0,1e10];
        [~,y] = ode15s(@h2o2_model,tspan,x0,options,N,GSH,GSSG,lambda,...
            PA,k_met,Vc,Ve,...
            k_kat,Km_kat,Ki_h2o2,alpha_max_kat(k),...
            k_ahp,Km_ahp,alpha_max_ahp,...
            alpha_max_grx,...
            n_oxyr,Km_oxyr,kon_oxyr_h2o2,kon_oxyr_grx,koff_oxyr_grx);
        Hin_rv(k,i) = y(end,2);
    end
end

%% Identify bistable region
tol = 1e-6;
bistability = zeros(length(alpha_max_kat),length(Hout));
for k=1:length(alpha_max_kat)
    for i=1:length(Hout)
        if (abs(Hin_fw(k,i)-Hin_rv(k,i))>tol)
            bistability(k,i)=1;
        end
    end
end

%% plot

figure();
hold on;

pcolor(bistability);

%% plot

figure();
hold on;

for k=1:length(alpha_max_kat)
    plot(Hout,Hin_fw(k,:),Hout,Hin_rv(k,:));
end
set(gca,'XScale','log');
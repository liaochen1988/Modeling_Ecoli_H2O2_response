clear all;
clc;

% HondropMatthews2004PlosBiology
time_hondrop2004 = [0,3,5,7,10,13,16,21,30,45,60,90]; % min
metE_act_hondrop2004 = [100.00,83.12,72.10,61.63,51.01,39.02,32.39,23.27,12.61,6.37,4.51,3.22]/100;

%% optimization
options = optimoptions('lsqnonlin','Display','iter','TolX',1e-6,'TolFun',1e-6,'MaxFunEval',Inf);
p0 = [0.000000001]; % kon
lb = [0];
ub = [Inf];
[sol,~,~,exitflag] = lsqnonlin(@mse_metE_decay,p0,lb,ub,options,time_hondrop2004,metE_act_hondrop2004);
assert(exitflag>0);

kon = sol(1); % 1/uM/min 
koff = kon/1.4; % 1/uM/min

%% plot comparison between experiment and prediction
[~,tpred,ypred] = mse_metE_decay(kon,time_hondrop2004,metE_act_hondrop2004);

figure();
hold on;
plot(tpred,ypred(:,1)/ypred(1,1),'k-');
plot(time_hondrop2004,metE_act_hondrop2004,'ko');
axis square;
box on;
xlabel('Time (min)');
ylabel('MetE activity');
axis([0,100,-0.1,1.1]);

function [dy,tpred,ypred] = mse_metE_decay(p,xdata,ydata)

% K_mix = 1.4 in HondropMatthews2004PlosBiology
kon = p(1);
koff = kon/1.4;

GSSG = 5*1e3; % uM
MetE_SH = 50; % uM

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,4),'NonNegative',1:4);
tspan = [0,100]; % min
x0 = [MetE_SH,GSSG,0,0];
[tpred,ypred] = ode15s(@metE_glutathionylation_model,tspan,x0,options,kon,koff);

dy = pchip(tpred,ypred(:,1),xdata)/MetE_SH-ydata;

end
function dy = mse_vary_h2o2_foxyrox_dyn(p,xdata,ydata,n_oxyr,Km_oxyr,k_oxyr_ox,Keq)

k_oxyr_ox1 = 10^p(1);
k_oxyr_red1 = k_oxyr_ox1*Keq;

% fixed GSH and GSSG concentration
GSSG = 1e-4; % M
GSH = 0.025; % M
oxyr_total = 0.01; % uM

Hin = xdata;
ypred = zeros(1,length(xdata));
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,3),'NonNegative',1:3);

for k=1:length(Hin)
    % before shift
    x0 = [0,0,oxyr_total];
    tspan = [0,1e6];
    [~,x_bs] = ode15s(@oxyr_kinetic_model_in_vitro,tspan,x0,options,...
        k_oxyr_ox,n_oxyr,Km_oxyr,k_oxyr_ox1,k_oxyr_red1,GSH,GSSG);
    
    % after shift
    x0 = x_bs(end,:);
    x0(1) = Hin(k); % uM
    tspan = [0,30];
    [~,x_as] = ode15s(@oxyr_kinetic_model_in_vitro,tspan,x0,options,...
        k_oxyr_ox,n_oxyr,Km_oxyr,k_oxyr_ox1,k_oxyr_red1,GSH,GSSG);
    
    ypred(k) = x_as(end,2)/oxyr_total;
end

dy = [ypred-ydata]';

end
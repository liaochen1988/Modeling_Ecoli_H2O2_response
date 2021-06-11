function dy = mse_vary_h2o2_foxyrox_ss(p,xdata1,ydata1,n_oxyr,Km_oxyr,k_oxyr_ox,Keq)

k_oxyr_ox1 = 10^p(1);
k_oxyr_red1 = k_oxyr_ox1*Keq;

% fixed GSH and GSSG concentration
GSSG = 1e-4; % M
GSH = 0.025; % M
Grx = 10*1e-6; % M
Hin = xdata1;
kswitch = k_oxyr_ox*Hin.^n_oxyr./(Km_oxyr^n_oxyr+Hin.^n_oxyr);
ypred1 = (kswitch+k_oxyr_ox1*Grx*GSSG.^4)./(kswitch+k_oxyr_red1*Grx*GSH.^8+k_oxyr_ox1*Grx*GSSG.^4);

dy = [ypred1-ydata1]';

end
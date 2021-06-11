function dy = mseHill_v4(p,xdata1,ydata1,hill_opt_1,Keq)

k_oxyr_ox1 = 10^p(1);
k_oxyr_red1 = k_oxyr_ox1*Keq;

n_oxyr = hill_opt_1(1);
Km_oxyr = hill_opt_1(2);
k_oxyr_ox = hill_opt_1(3);

% fixed GSH and GSSG concentration
GSSG = 0.1; % mM
GSH = 25; % mM
Hin = xdata1;
kswitch = k_oxyr_ox*Hin.^n_oxyr./(Km_oxyr^n_oxyr+Hin.^n_oxyr);
ypred1 = (kswitch+k_oxyr_ox1*GSSG.^4)./(k_oxyr_red1*GSH.^8+kswitch+k_oxyr_ox1*GSSG.^4);

dy = [ypred1-ydata1]';

end
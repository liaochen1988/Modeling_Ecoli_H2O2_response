function dy = mse_vary_gshgssgr_foxyrox(p,xdata,ydata)

% varied GSH and GSSG concentration
GSHGSSGR = xdata;
ypred1 = 1./(1+10^p*GSHGSSGR.^4);

dy = [ypred1-ydata]';

end
function dy = mseHill_v5(p,xdata1,ydata1)

% varied GSH and GSSG concentration
GSHGSSGR = xdata1;
ypred1 = 1./(1+10^p*GSHGSSGR.^4);

dy = [ypred1-ydata1]';

end
function dy = mseHill_v3(p,xdata,ydata,hill_opt_1)
ypred = xdata.^hill_opt_1(1)./(xdata.^hill_opt_1(1)*(1+p(1)/hill_opt_1(3))+p(1)/hill_opt_1(3)*hill_opt_1(2)^hill_opt_1(1));
dy = ypred-ydata;
end
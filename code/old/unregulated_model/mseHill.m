function dy = mseHill(p,xdata,ydata)
ypred = xdata.^p(1)./(xdata.^p(1)+p(2)^p(1));
dy = ypred-ydata;
end
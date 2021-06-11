function dy = mse_vary_h2o2_kobs(p,xdata,ydata)

ypred = p(3)*xdata.^p(1)./(xdata.^p(1)+p(2)^p(1));
dy = ypred-ydata;

end
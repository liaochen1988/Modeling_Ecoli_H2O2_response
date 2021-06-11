function dxdt = oxyr_kinetic_model_in_vitro(t,x,kon_oxyr_h2o2,n_oxyr,Km_oxyr,kon_oxyr_grx,koff_oxyr_grx,GSH,GSSG,Grx)

h2o2 = x(1);
oxyr_ox = x(2);
oxyr_red = x(3);

kswitch = kon_oxyr_h2o2*h2o2^n_oxyr/(Km_oxyr^n_oxyr+h2o2^n_oxyr);
dxdt(1) = -n_oxyr*kswitch*oxyr_red;
dxdt(2) = kswitch*oxyr_red + (kon_oxyr_grx*GSSG^4*oxyr_red - koff_oxyr_grx*GSH^8*oxyr_ox)*Grx;
dxdt(3) = -dxdt(2);

dxdt = real(dxdt');

end


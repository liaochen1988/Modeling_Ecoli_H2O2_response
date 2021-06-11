function dxdt = oxyr_kinetic_model_in_vitro(t,x,k_oxyr_ox,n_oxyr,Km_oxyr,k_oxyr_ox1,k_oxyr_red1,GSH,GSSG)

h2o2 = x(1);
oxyr_ox = x(2);
oxyr_red = x(3);

kswitch = k_oxyr_ox*h2o2^n_oxyr/(Km_oxyr^n_oxyr+h2o2^n_oxyr);
dxdt(1) = -n_oxyr*kswitch*oxyr_red;
dxdt(2) = kswitch*oxyr_red + k_oxyr_ox1*GSSG^4*oxyr_red - k_oxyr_red1*GSH^8*oxyr_ox;
dxdt(3) = -dxdt(2);

dxdt = real(dxdt');

end


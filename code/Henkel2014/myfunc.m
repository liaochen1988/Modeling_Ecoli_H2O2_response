function dxdt = myfunc(t,x,D,cin_glc,vmax_glc,Km_glc,vin_o2_100)

cx = x(1);
c_glc = x(2);
c_o2 = x(3);
c_qh2 = x(4);
c_q = x(5);
c_e2h2 = x(6);
a_oxi = x(7);
a_dh = x(8);

v_glc = vmax_glc*c_glc/(Km_glc+c_glc);
a = vin_o2/vin_o2_100;

dxdt = zeros(6,1);

dxdt(1) = (mu-D)*cx;
dxdt(2) = D*cin_glc-v_glc*cx-D*c_glc;
dxdt(3) = vin_o2-0.5*v_oxi*cx-D*c_o2;
dxdt(4) = -v_oxi+v_dh+vsyn_qh2-mu*c_qh2;
dxdt(5) = v_oxi-v_dh-mu*c_q;
dxdt(6) = -v_dh+12*v_glc-v_mu-v_ferm;
dxdt(7) = v_syn_oxi-mu*a_oxi;
dxdt(8) = v_syn_dh-mu*a_dh;

end


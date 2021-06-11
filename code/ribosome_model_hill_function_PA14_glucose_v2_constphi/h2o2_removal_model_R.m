function [dxdt,r_d,r_g] = h2o2_removal_model_R(t,x,phi,k_o,k_u,k_d,Km_o,Km_u,Km_d,Ki_u,f_o,f_u,n_u,e_g,k_diff,Vin,Vout,r_g0,r_mc,...
                                               beta,m_r,m_e,k_r,K_ma,K_ix,k_x,d_x,K_ia)

% variables
h2o2_in = x(1);
h2o2_out = x(2);
carbon = x(3);
ncell = x(4);
aa = x(5);
Rprot = x(6);
Xreg = x(7);
Eprot = x(8);

% rate equations
%r_o = k_o*carbon/(Km_o+carbon); % carbon oxidation rate
r_o = k_o*carbon;
r_a = k_r*Rprot*aa/(K_ma+aa);
activity = 1/(1+(h2o2_in/Ki_u)^n_u);
r_u = k_u*Eprot*phi*activity*carbon/(Km_u+carbon)*K_ia/(K_ia+aa); % carbon uptake rate, producing energy
r_d = k_d*Eprot*(1-phi)*activity*h2o2_in/(Km_d+h2o2_in)*carbon/(Km_u+carbon); % H2O2 consumption rate, consuming energy
r_g = r_a/beta; % bacterial growth rate
r_r = r_a/m_r*K_ix/(K_ix+Xreg);
r_e = r_a/m_e*(1-K_ix/(K_ix+Xreg));
r_x = k_x*K_ma/(K_ma+aa);

% differential equations
dxdt = zeros(7,1);
dxdt(1) = f_u*(r_g+r_g0)-r_d-k_diff*(h2o2_in-h2o2_out);
dxdt(2) = f_o*r_o+ncell*k_diff*(h2o2_in-h2o2_out)*Vin/(Vout-ncell*Vin);
dxdt(3) = -r_o-r_u*ncell*Vin/Vout;
dxdt(4) = (r_g-r_mc)*ncell;
dxdt(5) = r_u*e_g*beta-r_a-r_g*aa;
dxdt(6) = r_r-r_g*Rprot;
dxdt(7) = r_x-d_x*Xreg-r_g*Xreg;
dxdt(8) = r_e-r_g*Eprot;

end


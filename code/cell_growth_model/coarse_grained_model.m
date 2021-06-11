function dx = coarse_grained_model(x,k_n,k_f,k_t,rho,alpha_b,A_etc,k_atp,ca,cb,etc_max,nf)

d_phi_r = x(1);
d_phi_f = x(2);
d_phi_b = x(3);
lambda = x(4);
f_etc = x(5);

% fluxes
J_b = d_phi_b*nu_b; % flux prossed by ribosomes/biomass
J_r = d_phi_r*nu_r;
J_f = d_phi_f*nu_f;


% surface area/volume
SVr = ca-cb*lambda;

dx = zeros(5,1);
dx(1) = k_n*d_phi_r-k_f*d_phi_f-k_t*d_phi_b;
dx(2) = 1-d_phi_r-d_phi_f-d_phi_b;
dx(3) = lambda-k_t*d_phi_b/rho;

if (etc_max < alpha_b*lambda*f_etc*A_etc/k_atp/SVr)
    dx(4) = etc_max-alpha_b*lambda*f_etc*A_etc/k_atp/SVr;
    dx(5) = nf*k_f*d_phi_f-alpha_b*lambda*(1-f_etc);
else
    dx(4) = 0;
    dx(5) = nf*k_f*d_phi_f-0;
end

end


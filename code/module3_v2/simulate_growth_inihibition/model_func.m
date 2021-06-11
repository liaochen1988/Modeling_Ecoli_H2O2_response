function dy = model_func(t,x,N,k_met,PA,k0_cat,k0_ahp,Km_ahp,Vc,Ve,alpha_cat,alpha_ahp,lambda,n_oxyr,k_oxyr_red,k_oxyr_ox,Km_oxyr)

Ho = x(1);
Hi = x(2);
Cat = x(3);
Ahp = x(4);

% fraction of oxidized OxyR 
f = Hi^n_oxyr/(Hi^n_oxyr*(1+k_oxyr_red/k_oxyr_ox)+k_oxyr_red/k_oxyr_ox*Km_oxyr^n_oxyr);

dy = zeros(4,1);
dy(1) = -(Ho-Hi)*PA/Ve*N;
dy(2) = (k_met+(Ho-Hi)*PA)/Vc-k0_cat*Cat*Hi-k0_ahp*Ahp*Hi/(Km_ahp+Hi);
dy(3) = alpha_cat*f/Vc-Cat*lambda;
dy(4) = alpha_ahp*f/Vc-Ahp*lambda;

end
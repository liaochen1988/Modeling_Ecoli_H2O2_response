function dy = model_func_v2(t,x,Ho,k_met,PA,k0_cat,k0_ahp,Km_ahp,Vc,Ve,alpha_cat,alpha_ahp,lambda,n_oxyr,k_oxyr_red,k_oxyr_ox,Km_oxyr)

Hi = x(1);
Cat = x(2);
Ahp = x(3);

% fraction of oxidized OxyR 
f = Hi^n_oxyr/(Hi^n_oxyr*(1+k_oxyr_red/k_oxyr_ox)+k_oxyr_red/k_oxyr_ox*Km_oxyr^n_oxyr);

dy = zeros(3,1);
dy(1) = (k_met+(Ho-Hi)*PA)/Vc-k0_cat*Cat*Hi-k0_ahp*Ahp*Hi/(Km_ahp+Hi);
dy(2) = alpha_cat*f/Vc-Cat*lambda;
dy(3) = alpha_ahp*f/Vc-Ahp*lambda;

end
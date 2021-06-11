function dy = h2o2_model(t,x,N,GSH,GSSG,lambda,...
                         PA,k_met,Vc,Ve,...
                         k_cat,Km_kat,Ki_h2o2,alpha_max_kat,...
                         k_ahp,Km_ahp,alpha_max_ahp,...
                         alpha_max_grx,...
                         n_oxyr,Km_oxyr,kon_oxyr_h2o2,kon_oxyr_grx,koff_oxyr_grx)

Ho = x(1);
Hi = x(2);
Kat = x(3);
Ahp = x(4);
Grx = x(5);

% fraction of oxidized OxyR 
kswitch = kon_oxyr_h2o2*Hi.^n_oxyr./(Km_oxyr^n_oxyr+Hi.^n_oxyr);
f = (kswitch+kon_oxyr_grx*Grx*GSSG.^4)./(kswitch+koff_oxyr_grx*Grx*GSH.^8+kon_oxyr_grx*Grx*GSSG.^4);

dy = zeros(5,1);
dy(1) = -(Ho-Hi)*PA/Ve*N;
dy(2) = (k_met+(Ho-Hi)*PA)/Vc-k_cat*Kat*Hi/(Km_kat+Hi)/(1+Hi/Ki_h2o2)-k_ahp*Ahp*Hi/(Km_ahp+Hi);
dy(3) = alpha_max_kat*f-Kat*lambda;
dy(4) = alpha_max_ahp*f-Ahp*lambda;
dy(5) = alpha_max_grx*f-Grx*lambda;

end
function dy = mseOxyR_v3(p,k_met,PA,k_cat,k_ahp,n_ahp,hill_opt_invitro,h2o2_invivo,oxyr_ox_invivo)

options = optimset('Display','off','TolX',1e-12); % show iterations
Hout = 10.^[-9:0.05:-4]; % M
Hin = zeros(size(Hout));
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin(i+1);
    end
    
    [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,k_ahp,10^p(1),n_ahp);
    while Hin(i)<0
        x0 = 10.^(-9+rand*5);
        [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,k_ahp,10^p(1),n_ahp);
    end
    assert(exitflag>0);
end

n_oxyr = hill_opt_invitro(1);
Km_oxyr =  hill_opt_invitro(2);
kox_oxyr = hill_opt_invitro(3);
kred_oxyr = hill_opt_invitro(4);
oxyr_ox_invivo_pred = Hin.^n_oxyr./(Hin.^n_oxyr*(1+kred_oxyr/kox_oxyr)+kred_oxyr/kox_oxyr*Km_oxyr^n_oxyr);

dy = pchip(Hout,oxyr_ox_invivo_pred,h2o2_invivo)-oxyr_ox_invivo;

end


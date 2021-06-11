function dy = mseOxyR_kahp_P(p,k_met,Km_ahp,k_cat,n_ahp,hill_opt_invitro,h2o2_invivo,oxyr_ox_invivo)

options = optimset('Display','off','TolX',1e-12); % show iterations
Hout = 10.^[-9:0.05:-4]; % M
Hin = zeros(size(Hout));
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin(i+1);
    end
    
    [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),10^p(2),k_cat,10^p(1),Km_ahp,n_ahp);
    while Hin(i)<0
        x0 = 10.^(-9+rand*5);
        [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),10^p(2),k_cat,10^p(1),Km_ahp,n_ahp);
    end
    assert(exitflag>0);
end
oxyr_ox_invivo_pred = Hin.^hill_opt_invitro(1)./(Hin.^hill_opt_invitro(1)+hill_opt_invitro(2)^hill_opt_invitro(1));

dy = pchip(Hout,oxyr_ox_invivo_pred,h2o2_invivo)-oxyr_ox_invivo;

end


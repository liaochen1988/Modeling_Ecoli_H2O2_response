function dy = mseOxyR_v4(p,k_met,PA,k_cat,n_ahp,hill_opt_invitro,h2o2_invivo,oxyr_ox_invivo,h2o2_ahp_invivo,Jahp_invivo)

%% fraction of oxidized OxyR
options = optimset('Display','off','TolX',1e-12); % show iterations
Hout = 10.^[-9:0.05:-4]; % M
Hin = zeros(size(Hout));
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin(i+1);
    end
    
    [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,10^p(1),10^p(2),n_ahp);
    while Hin(i)<0
        x0 = 10.^(-9+rand*5);
        [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,10^p(1),10^p(2),n_ahp);
    end
    assert(exitflag>0);
end

n_oxyr = hill_opt_invitro(1);
Km_oxyr =  hill_opt_invitro(2);
kox_oxyr = hill_opt_invitro(3);
kred_oxyr = hill_opt_invitro(4);
oxyr_ox_invivo_pred = Hin.^n_oxyr./(Hin.^n_oxyr*(1+kred_oxyr/kox_oxyr)+kred_oxyr/kox_oxyr*Km_oxyr^n_oxyr);

%% H2O2 decomposition rate
options = optimset('Display','off','TolX',1e-12); % show iterations
Hout = 10.^[-9:0.05:-4]; % M
Hin = zeros(size(Hout));
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin(i+1);
    end
    
    [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,0,Hout(i),PA,0,10^p(1),10^p(2),n_ahp);
    while Hin(i)<0
        x0 = 10.^(-9+rand*5);
        [Hin(i),~,exitflag] = fzero(@calcHin,x0,options,0,Hout(i),PA,0,10^p(1),10^p(2),n_ahp);
    end
    assert(exitflag>0);
end
Jahp = pchip(Hout,10^p(1).*Hin./(Hin+10^p(2)),h2o2_ahp_invivo);

%% calculate difference
%dy1 = (pchip(Hout,oxyr_ox_invivo_pred,h2o2_invivo)-oxyr_ox_invivo)/max(oxyr_ox_invivo);
dy1 = (pchip(Hout,oxyr_ox_invivo_pred,h2o2_invivo)-oxyr_ox_invivo);
dy2 = Jahp/max(Jahp)-Jahp_invivo/max(Jahp_invivo);
dy = [dy1,dy2]';

end


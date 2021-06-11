clear all;
clc;

% Use values from SeaverImlay2001JB for 
% P (permeability)
% A (area) 
% k_met (H2O2 production rate by metabolism)
% k_cat (H2O2 degradation rate of catalase)
P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20; % mol/s
k_cat = 1.05e-15; % L/s
Km_kat = 3.9e-3; % M
Ki_katG_h2o2 = 19.55/1e6; % M, obtained from independent fitting
n_ahp = 1;

n_oxyr = 1.3614;
Km_oxyr = 41.3179/1e6; % M
kon_oxyr_h2o2 = 9.9863;
kon_oxyr_grx = 1.4910e7;
koff_oxyr_grx = 1.5218e15;

GSSG = 1e-4; % M
GSH = 0.025; % M
Grx = 10*1e-6; % M

load('par_est_ahp_kinetics.mat');

%% change Km_ahp

figure();
hold on;

Km_ahp_arr = 10.^[-9:1:-4];
Hout = 10.^[-9:0.001:-4]; % M
cc = jet(length(Km_ahp_arr));
options = optimset('Display','off','TolX',1e-12); % show iterations

for k=length(Km_ahp_arr):-1:1
    Hin_f = zeros(size(Hout));
    for i=length(Hout):-1:1
        if(i==length(Hout))
            x0 = Hout(i);
        else
            x0 = Hin_f(i+1);
        end
        [Hin_f(i),fval,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,Km_kat,Ki_katG_h2o2,k_ahp,Km_ahp_arr(k),n_ahp);
        while Hin_f(i)<0
            x0 = 10.^(-9+rand*5);
            [Hin_f(i),~,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,Km_kat,Ki_katG_h2o2,k_ahp,Km_ahp_arr(k),n_ahp);
        end
        assert(exitflag>0);
    end
    
    kswitch =  kon_oxyr_h2o2*Hin_f.^n_oxyr./(Km_oxyr^n_oxyr+Hin_f.^n_oxyr);
    frac_oxyr_ox_pred = (kswitch+kon_oxyr_grx*Grx*GSSG^4)./(kswitch+kon_oxyr_grx*Grx*GSSG^4+koff_oxyr_grx*Grx*GSH^8);
    plot(Hout*1e6,frac_oxyr_ox_pred,'k-','Color',cc(k,:));
end

axis([1e-3,1e2,-0.1,1.1]);
axis square;
box on;
ylabel('Oxidized OxyR (%)');
xlabel('Extracellular H_2O_2 (\muM)');
set(gca,'XScale','log');
set(gca,'XTick',10.^[-3:1:2]);
set(gca,'YTick',[0:0.2:1.0]);
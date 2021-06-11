clear all;
clc;

%% parameters
P=1.6e-3; % cm/s
A=1.41e-7; % cm^2
PA = P*A/1e3; % L/s
k_met = 4.5e-20; % mol/s
%k_met = 0;
k_cat = 2.7e-13; % L/s
k_ahp = 2.1e-18; % mol/s
Km_ahp = 1.2e-10; % M
n_ahp = 1;
Km_oxyr = 8.8342e-08;

%% Predict OxyR oxidation percentage
Hout = 10.^[-12:0.01:-4]; % M
Hin = zeros(size(Hout));
options = optimset('Display','off','TolX',1e-12); % show iterations
for i=length(Hout):-1:1
    if(i==length(Hout))
        x0 = Hout(i);
    else
        x0 = Hin(i+1);
    end
    [Hin(i),fval,exitflag] = fzero(@calcHin,x0,options,k_met,Hout(i),PA,k_cat,k_ahp,Km_ahp,n_ahp);
    assert(exitflag>0);
end

% simulated derivatives
pp=spline(Hout,Hin);
p_der=fnder(pp,1);
deriv = ppval(p_der,Hout).*Hout;
deriv = deriv*Km_oxyr./(Km_oxyr+Hin).^2;

% analytic derivatives
Delta = Km_ahp*PA+Km_ahp*k_cat-Hout*PA-k_met+k_ahp;
Gamma = Km_ahp*PA+Km_ahp*k_cat+Hout*PA+k_met-k_ahp;
Eta = 4*Km_ahp*(PA+k_cat)*(k_met+PA*Hout);
%deriv_ana = Hout.*(1+Gamma./sqrt(Delta.^2+Eta))*PA/2/(PA+k_cat);
%deriv_ana = Hout.*(1+Gamma./(Delta+Eta./Delta/2))*PA/2/(PA+k_cat);
deriv_ana = Km_ahp*PA*Hout*(Km_ahp*(PA+k_cat)+k_ahp)./(Delta.^2+Eta/2);
deriv_ana = deriv_ana*Km_oxyr./(Km_oxyr+Hin).^2;

%% plot derivatives
figure();
hold on;

plot(Hout*1e6, deriv, 'r-');
plot(Hout*1e6, deriv_ana,'b-');
%legend('In vivo','In vitro');
%title('Simulation');
xlim([1e0,1e2]);
%axis([1e-3,1e2,-0.1,1.1]);
axis square;
box on;
% ylabel('Oxidized OxyR (%)');
% xlabel('Extracellular H_2O_2 (\muM)');
set(gca,'XScale','log');
% set(gca,'XTick',10.^[-3:1:2]);
% set(gca,'YTick',[0:0.2:1.0]);

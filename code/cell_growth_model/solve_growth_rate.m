clear all;
clc;

k_f = 4;
k_t = 4;
rho = 0.76;
alpha_b = 0.5;
A_etc = 1.0;
k_atp = 270; % ATP?s
ca = 9.3; 
cb = 2.8;
etc_max = 1;
nf = 2;

x0 = rand(5,1);
options = optimoptions('fsolve','Display','off','TolX',1e-6);

k_n = [0:0.1:1];
sol = zeros(length(k_n),5);
for i=1:length(k_n)
    i
    [sol(i,:),~,exitflag] = fsolve(@coarse_grained_model,x0,options,k_n(i),...
                                    k_f,k_t,rho,alpha_b,A_etc,k_atp,ca,cb,etc_max,nf);
    assert(exitflag>0);
end

figure();
hold on;

plot(sol(:,4),sol(:,1),'r-');
plot(sol(:,4),sol(:,2),'b-');
plot(sol(:,4),sol(:,3),'g-');

axis square;
box on;
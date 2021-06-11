clear all;
clc;

k_oxyr_ox = 9.9863; % 1/s
n_oxyr = 1.3614; 
Km_oxyr = 41.3179; % uM
k_oxyr_ox1 = 149.1049; % 1/s/M^4
k_oxyr_red1 = 1.5218e+10; % 1/s/M^8
oxyr_total = 0.01; % uM
GSH = 0.025; % M
GSSG = 1e-4; % M

options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,3),'NonNegative',1:3);

% before shift
x0 = [0,0,oxyr_total];
tspan = [0,1e10];
[t_bs,x_bs] = ode15s(@oxyr_kinetic_model_in_vitro,tspan,x0,options,...
    k_oxyr_ox,n_oxyr,Km_oxyr,k_oxyr_ox1,k_oxyr_red1,GSH,GSSG);

% after shift
x0 = x_bs(end,:);
x0(1) = 0.2; % uM
tspan = [0,10000];
[t_as,x_as] = ode15s(@oxyr_kinetic_model_in_vitro,tspan,x0,options,...
    k_oxyr_ox,n_oxyr,Km_oxyr,k_oxyr_ox1,k_oxyr_red1,GSH,GSSG);

%figure();

subplot(2,1,1);
hold on;
plot([t_bs-max(t_bs);t_as]/60,[x_bs(:,1);x_as(:,1)]);
box on;
xlabel('Time (min)');
ylabel('[H_2O_2] (\muM)');
axis([-1,60,0,2]);

subplot(2,1,2);
hold on;
plot([t_bs-max(t_bs);t_as]/60,[x_bs(:,2);x_as(:,2)]/oxyr_total);
box on;
xlabel('Time (min)');
ylabel('Oxidized OxyR fraction');
axis([-1,60,-0.1,1.1]);

kswitch = k_oxyr_ox*x0(1)^n_oxyr/(Km_oxyr^n_oxyr+x0(1)^n_oxyr);
foxyrox_ss = (kswitch+k_oxyr_ox1*GSSG^4)/(k_oxyr_red1*GSH^8+kswitch+k_oxyr_ox1*GSSG^4);
plot([-1,60],[foxyrox_ss,foxyrox_ss],'k--');
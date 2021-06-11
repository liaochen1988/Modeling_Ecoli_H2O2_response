clear all;
clc;

Gor = 2.8499 * 1e-6; % M
NADPH = 1.23e-4;     % M
GSSG0 = 0.00237;     % M
fc = [1,2,5,10,20,50,100];
cls = jet(length(fc));
options=odeset('RelTol',1e-10,'AbsTol',1e-10, 'NonNegative',1);

figure();
hold on;

for i=1:length(fc)
    [t,x] = ode15s(@my_func,[0,60],NADPH,options,GSSG0 * fc(i),Gor);
    
    plot(t,x(:,1)/x(1,1),'k-','Color',cls(i,:));
end

axis([-1,60,-0.1,1.5]);
xlabel('Time (s)');
ylabel('Relative NADPH');

t_1mM = [0, 5, 15, 30, 60]; % sec
NADPH_1mM = [1.000	1.163	0.645	0.697	0.568]; % Relative
plot(t_1mM, NADPH_1mM, 'ko');
axis square;
box on;
clear all;
clc;

hp = readParameters();

%% Simulate different H2O2 concentration

kcat = 0.12;        % s
GSH  = 16583;       % uM
GSSG = 2400;        % uM
Grx  = 10;          % uM

H2O2 = 10.^[0:0.1:9]; % uM
GrowthRate = zeros(size(H2O2));
Met = zeros(size(H2O2));
Metox = zeros(size(H2O2));

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3,4,5,6],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_PMC);

%   time span
tspan = [0 10^10];

%   initial condition
x0 = [10,0,0,10,0,10];
tic;
[~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,kcat,0,GSH,GSSG,Grx,hp);
if (isempty(te))
    x0 = x(end,:);
else
    error('Error: Oscillation Detected for the current parameter set!');
end

for i=1:length(H2O2)
    tic;
    [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,kcat,H2O2(i),GSH,GSSG,Grx,hp);
    if (~isempty(te))
        error('Error: Oscillation Detected for the current parameter set!');
    end
    x = x(end,:);   %   only keep the steady state solution
    x0 = x;
    
    Met(i)   = x0(1);
    Metox(i)       = x0(2);
    [~,GrowthRate(i)] = Ecoli_GR_ODE_PMC(0,x0,kcat,H2O2(i),GSH,GSSG,Grx,hp);
end

%figure();

subplot(1,2,1);
hold on;
plot(H2O2/1e3, GrowthRate/GrowthRate(1));
set(gca,'XScale','log');
xlabel('H2O2 (mM)');
ylabel('Growth rate (h^{-1})');
axis square;
box on;
axis([1e-3,1e6,-0.1,1.1]);

subplot(1,2,2);
hold on;
plot(H2O2/1e3, Metox./(Metox + Met));
set(gca,'XScale','log');
xlabel('H2O2 (mM)');
ylabel('Fraction of oxidized Met');
axis square;
box on;
axis([1e-3,1e6,-0.1,1.1]);


clear all;
clc;

hp = readParameters();

%% Simulate different GSH/GSSG ratio

kcat = 0.12;        % s
GSH  = [1,5,10,15,20] * 1e3;       % uM
Ratio = 10.^[-2:0.1:2];
Grx  = 10.^[-2:1:2];       % uM

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3,4,5,6],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_PMC);

%   time span
tspan = [0 10^10];

figure();
for k=1:length(Grx)
    for i=1:length(GSH)
        
        GrowthRate = zeros(size(Ratio));
        MetE = zeros(size(Ratio));
        MetEox = zeros(size(Ratio));
        
        x0 = [10,0,0,10,0,10]; %   initial condition
        for j=1:length(Ratio)
            [k,i,j]
            GSSG = GSH(i)/Ratio(j);
            tic;
            [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,kcat,0,GSH(i),GSSG,Grx(k),hp);
            if (~isempty(te))
                error('Error: Oscillation Detected for the current parameter set!');
            end
            x = x(end,:);   %   only keep the steady state solution
            x0 = x;
            
            MetE(j)   = x0(4);
            MetEox(j) = x0(5);
            [~,GrowthRate(j)] = Ecoli_GR_ODE_PMC(0,x0,kcat,0,GSH(i),GSSG,Grx(k),hp);
        end
        
        subplot(2,5,k);
        hold on;
        plot(Ratio, GrowthRate * 3600);
        set(gca,'XScale','log');
        xlabel('GSH/GSSG');
        ylabel('Growth rate (h^{-1})');
        axis square;
        box on;
        %axis([1e-3,1e6,-0.1,1.1]);
        
        subplot(2,5,k+5);
        hold on;
        plot(Ratio, MetEox./(MetEox + MetE));
        set(gca,'XScale','log');
        xlabel('GSH/GSSG');
        ylabel('Fraction of oxidized MetE');
        axis square;
        box on;
        %axis([1e-3,1e6,-0.1,1.1]);
    end
end

%% Simulate different GSH/GSSG ratio

kcat = 0.12;        % s
GSH  = 5 * 1e3;     % uM
Ratio = 10.^[-3:0.1:1];
Grx  = 10.^[-2:1:2];       % uM

%   set Matlab solvers
tol = 1e-6;
options_ode15s_PMC = odeset('NonNegative',[1,2,3,4,5,6],...
    'RelTol',tol,...
    'AbsTol',tol,...
    'Events',@myEvent_PMC);

%   time span
tspan = [0 10^10];

figure();
for k=1:length(Grx)
    GrowthRate = zeros(size(Ratio));
    MetE = zeros(size(Ratio));
    MetEox = zeros(size(Ratio));
    
    x0 = [10,0,0,10,0,10]; %   initial condition
    for j=length(Ratio):-1:1
        [k,i,j]
        GSSG = GSH/Ratio(j);
        tic;
        [~,x,te] = ode15s(@Ecoli_GR_ODE_PMC,tspan,x0,options_ode15s_PMC,kcat,0,GSH,GSSG,Grx(k),hp);
        if (~isempty(te))
            error('Error: Oscillation Detected for the current parameter set!');
        end
        x = x(end,:);   %   only keep the steady state solution
        x0 = x;
        
        MetE(j)   = x0(4);
        MetEox(j) = x0(5);
        [~,GrowthRate(j)] = Ecoli_GR_ODE_PMC(0,x0,kcat,0,GSH,GSSG,Grx(k),hp);
    end
    
    subplot(1,2,1);
    hold on;
    plot(Ratio, GrowthRate * 3600);
    set(gca,'XScale','log');
    xlabel('GSH/GSSG');
    ylabel('Growth rate (h^{-1})');
    axis square;
    box on;
    %axis([1e-3,1e6,-0.1,1.1]);
    
    subplot(1,2,2);
    hold on;
    plot(Ratio, MetE./(MetEox + MetE));
    set(gca,'XScale','log');
    xlabel('GSH/GSSG');
    ylabel('MetE activity');
    axis square;
    box on;
    %axis([1e-3,1e6,-0.1,1.1]);
end
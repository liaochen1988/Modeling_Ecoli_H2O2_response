clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose 10 g/L
GLC_ex_input = 10.^[-5:0.1:0];

%   scaled dissolved oxygen level
a1 = (0.01:0.04:0.98);
a2 = (1.0:1.0:15.0);
a3 = (20.0:10.0:50.0);
a = [a1 a2 a3];

Num_of_State_Variable = 36;
Num_of_Flux = 11;
options=odeset('RelTol',1e-6,'AbsTol',1e-20, 'NonNegative',[1:Num_of_State_Variable]);

growth_rate = zeros(length(GLC_ex_input),length(a));

parfor i=1:length(GLC_ex_input)
    IC = Initial_Concentration(GLC_ex_input(i));
    
    gr_temp = zeros(1,length(a));
    for j=1:length(a)
        
        j
        
        called_by_ode15s = true;
        [t,x] = ode15s(@Kinetic_model,[-10.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a(j),called_by_ode15s);
        
        called_by_ode15s = false;
        gr_temp(j) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a(j),called_by_ode15s);
        
        IC = x(end,:);
    end
    growth_rate(i,:) = gr_temp*3600;    %   convert 1/s to 1/h
    
    plot(a,gr_temp);
    drawnow
    hold on;
end

axis square;

%%  plot
figure();

subplot(1,3,1);
cls = jet(length(a));
for j=1:length(a)
   plot(GLC_ex_input,growth_rate(:,j),'Color',cls(j,:));
   hold on;
end
axis square;
box on;
set(gca,'XScale','log');

subplot(1,3,2);
cls = jet(length(GLC_ex_input));
for j=1:length(GLC_ex_input)
   plot(a,growth_rate(j,:),'Color',cls(j,:));
   hold on;
end
axis square;
box on;

subplot(1,3,3);
[xdata,ydata]=meshgrid(GLC_ex_input,a);
pcolor(xdata,ydata,growth_rate');
shading interp
axis square;
box on;

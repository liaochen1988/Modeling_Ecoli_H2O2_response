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

growth_rate = zeros(length(a),length(GLC_ex_input));
parfor i=1:length(a)
    
    gr_temp = zeros(1,length(GLC_ex_input));
    IC = Initial_Concentration(GLC_ex_input(1));
    for j=1:length(GLC_ex_input)
        j
        IC(2) = GLC_ex_input(j);
        called_by_ode15s = true;
        [t,x] = ode15s(@Kinetic_model,[-10.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a(i),called_by_ode15s);
        
        called_by_ode15s = false;
        gr_temp(j) = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a(i),called_by_ode15s);
        
        IC = x(end,:);
    end
    growth_rate(i,:) = gr_temp*3600;    %   convert 1/s to 1/h
end


%%  plot

load('glucose-oxygen-growth-data-scanning-foreach-oxygen');
growth_rate_cp = growth_rate;

growth_rate_cp(growth_rate_cp<-1e-10) = 1e-6;

%   correcting negative growth rate values
for j=1:size(growth_rate,2)

    %   find the first place that a positive value is found

    if (growth_rate_cp(1,j)>1e-10)
        continue;
    end
    
    pos=0;
    for i=1:size(growth_rate,1)
        if (growth_rate_cp(i,j) > 1e-10)
            pos=i;
            break;
        end
    end
    
    %growth_rate_cp(1:pos-1,j) = pchip(a(pos:end),growth_rate(pos:end,j),a(1:pos-1));
    growth_rate_cp(1:pos-1,j) = growth_rate_cp(pos,j);
end

figure();

subplot(1,3,1);
cls = jet(length(a));
for j=IDX:length(a)
   plot(GLC_ex_input,growth_rate_cp(j,:),'Color',cls(j,:));
   hold on;
end
axis square;
box on;
set(gca,'XScale','log');

subplot(1,3,2);
cls = jet(length(GLC_ex_input));
for j=1:length(GLC_ex_input)
   plot(a,growth_rate_cp(:,j),'Color',cls(j,:));
   hold on;
end
axis square;
box on;
set(gca,'XScale','log');

subplot(1,3,3);
[gx,gy]=meshgrid(a,GLC_ex_input);
surf(gx,gy,growth_rate_cp');
shading interp
axis square;
box on;
set(gca,'XScale','log');
set(gca,'YScale','log');

%%  fitting with a simple glucose-oxygen model
xdata1 = reshape(repmat(a',1,length(GLC_ex_input)),length(a)*length(GLC_ex_input),1)*0.214375/100*1000;   %   oxygen (uM)
xdata2 = reshape(repmat(GLC_ex_input,length(a),1),length(a)*length(GLC_ex_input),1)/180.156*1e+6;         %   glucose (uM)
xdata = [xdata1';xdata2']';
ydata = reshape(growth_rate_cp,length(a)*length(GLC_ex_input),1);

%   Hill coefficient 1
% F = @(p,xdata) (p(1)*(xdata(:,2)/p(2))+p(3).*(xdata(:,1)/p(4)).*(xdata(:,2)/p(2)))./(1+(xdata(:,1)/p(4))+(xdata(:,2)/p(2))+(xdata(:,1)/p(4)).*(xdata(:,2)/p(2)));
% p0 = [0.3,1e-4,0.9,2];
% [sol,resnorm,~,exitflag] = lsqcurvefit(F,p0,xdata,ydata);
% p=sol;
% ypred = (p(1)*(xdata(:,2)/p(2))+p(3).*(xdata(:,1)/p(4)).*(xdata(:,2)/p(2)))./(1+(xdata(:,1)/p(4))+(xdata(:,2)/p(2))+(xdata(:,1)/p(4)).*(xdata(:,2)/p(2)));
% ypred = reshape(ypred,length(a),length(GLC_ex_input));

%   Hill coefficient as a fitting parameter
F = @(p,xdata) (p(1)*(xdata(:,2)/p(2)).^p(5)+p(3).*(xdata(:,1)/p(4)).^p(6).*(xdata(:,2)/p(2)).^p(5))./(1+(xdata(:,1)/p(4)).^p(6)+(xdata(:,2)/p(2)).^p(5)+(xdata(:,1)/p(4)).^p(6).*(xdata(:,2)/p(2)).^p(5));
p0 = [0.3,0.55,0.9,4.28,1,1];
[sol,resnorm,~,exitflag] = lsqcurvefit(F,p0,xdata,ydata);
p=sol;
ypred = (p(1)*(xdata(:,2)/p(2)).^p(5)+p(3).*(xdata(:,1)/p(4)).^p(6).*(xdata(:,2)/p(2)).^p(5))./(1+(xdata(:,1)/p(4)).^p(6)+(xdata(:,2)/p(2)).^p(5)+(xdata(:,1)/p(4)).^p(6).*(xdata(:,2)/p(2)).^p(5));
ypred = reshape(ypred,length(a),length(GLC_ex_input));

[gx,gy]=meshgrid(a*0.214375/100*1000,GLC_ex_input/180.156*1e+6);
surf(gx',gy',ypred);
shading interp
hold on;
plot3(xdata1,xdata2,ydata,'ko','MarkerFaceColor','c');
axis square;
box on;
set(gca,'XScale','log');
set(gca,'YScale','log');
view([28,16]);
xlabel('Oxygen');
ylabel('Glucose');
zlabel('Growth rate');
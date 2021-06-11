clear all;

%   whether simulating fnr and arcA mutants
fnr_mutant = false;
arcA_mutant = false;

%   initial external glucose 10 g/L
GLC_ex_input = 10.^[-5:0.2:0];

%   scaled dissolved oxygen level
a = 40; % aerobic condition

MET           = zeros(length(GLC_ex_input),1);
RIB           = zeros(length(GLC_ex_input),1);
ppGpp         = zeros(length(GLC_ex_input),1);
growthRate    = zeros(length(GLC_ex_input),1);

%% integrated model
integrated_model = true;
Num_of_State_Variable = 40;
options=odeset('RelTol',1e-6,'AbsTol',1e-6, 'NonNegative',[1:Num_of_State_Variable]);

c_met_e = 170;

IC = Initial_Concentration();
for i=1:length(GLC_ex_input)
    [i]
    IC(2) = GLC_ex_input(i);
    called_by_ode15s = true;
    [t,x] = ode15s(@Kinetic_model,[-1000.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,c_met_e);
    
    called_by_ode15s = false;
    res = Kinetic_model(0,x(end,:),arcA_mutant,fnr_mutant,a,called_by_ode15s,integrated_model,c_met_e);
    MET(i)          = res(1);
    RIB(i)          = res(2);
    ppGpp(i)        = res(3);
    growthRate(i)   = res(4) * 3600;
    
    IC = x(end,:);
end

%%  plot
figure();

%   color code
ColorPalette = [255,109,109;255,219,107;72,255,167;61,191,255;198,73,255]/255;

%   maker size
MS = 8;

%   line width
LW = 1.5;

%----------------------------------------------------------------
%  R protein fraction (affliated protein excluded) vs Growth rate
%----------------------------------------------------------------

subplot(1,3,1);
hold on;

l_rib = 7336*1.6;
total_aa = 3.00E6;
plot(growthRate, RIB * l_rib / total_aa / 1.6, 'k-', 'LineWidth', LW);
axis square;
box on;
xlabel('Growth rate (hr^{-1})');
ylabel('R fraction');

%   Paper 1: Modulation of chemical composition and other parameters of the
%   cell by growth rate
%   Hans Bremer & Patrick Dennis (1996)
growthRate_BP1996 = [0.6,1.0,1.5,2.0,2.5]*log(2);
phiR_BP1996       = [9,11,13.5,18.0,21.6]/100;
plot(growthRate_BP1996,phiR_BP1996,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(1,:));

%   Paper 2: Interdependence of cell growth and gene expression: origins and
%   consequences
%   Matthew Scott et al. (2010)
expData = [...
    [0.4,0.57,0.71,1.00,1.31,1.58];...                    %   growth rate
    [0.177,0.230,0.224,0.287,0.414,0.466]*0.76/1.6;...    %   R-protein fraction
    [0.03,0.02,0.03,0.05,0.07,0.15];...                   %   s.d. of growth rate
    [0.006,0.014,0.029,0.009,0.058,0.033]...              %   s.d. of R-protein fraction
    ]';
plot(expData(:,1),expData(:,2),'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(2,:));
errorxy(expData,'ColXe',3,'ColYe',4,'Marker','o','EdgeColor','black','MarkSize',MS,'ColorEB','k','EdgeEB',MS);

%   Paper 3: Quantifying the benefit of a proteome reserve in fluctuating
%   environments
%   Metteo Mori et al. (2017)
growthRate_MM2017 = [0.06,0.12,0.15,0.18,0.19,0.25,0.28,0.31,0.36,0.43,0.60,0.60,0.54,0.70,0.77,0.78,0.89,0.91,0.89,0.96,1.02,2.46];
phiR_MM2017       = [0.13,0.14,0.14,0.13,0.15,0.16,0.16,0.17,0.18,0.18,0.21,0.20,0.22,0.24,0.27,0.28,0.26,0.29,0.29,0.30,0.30,0.66]*0.76/1.6;
plot(growthRate_MM2017,phiR_MM2017,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(3,:));

%   Paper 4: Reduction of translating ribosomes enables Escherichi coli to
%   maintain elongation rates during slow growth
%   Xiongfeng Dai et al. (2017)
growthRate_XD2017 = [1.797,1.568,1.268,1.165,1.096,0.975,0.947,0.990,0.918,0.745,0.717,0.705,0.691,0.675,0.551,0.485,0.456,0.415,0.407,0.339,0.284,0.229,0.195,0.129,0.035,0];
phiR_XD2017       = [0.228,0.210,0.173,0.157,0.146,0.140,0.130,0.127,0.121,0.111,0.102,0.101,0.103,0.103,0.090,0.088,0.082,0.082,0.082,0.072,0.070,0.064,0.061,0.058,0.046,0.042];
plot(growthRate_XD2017,phiR_XD2017,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(4,:));

xlim([0,2]);
ylim([0,0.4]);

%----------------------------------------------------------------
%   ppGpp vs Growth rate
%----------------------------------------------------------------

subplot(1,3,2);
hold on;

plot(growthRate,ppGpp,'k-','LineWidth',LW);
axis square;
box on;
xlabel('Growth rate (hr^{-1})');
ylabel('ppGpp (\muM)');

%   conversion factor between OD460 and um^3
%   Paper: Molecular crowding limits translation and cell growth by S.
%   Klumpp et al. (2013)
convFac = 0.32;

%   Paper 1: Modulation of chemical composition and other parameters of the
%   cell by growth rate
%   Hans Bremer & Patrick Dennis (1996)
growthRate_BP1996 = [0.6,1.0,1.5,2.0,2.5]*log(2);
ppGpp_BP1996      = [55,38,22,15,10]/convFac;
plot(growthRate_BP1996,ppGpp_BP1996,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(1,:));

%   Paper 2: Control of rRNA and tRNA syntheses in Escherichi coli by guanosine tetraphosphate
%   J Ryals, R Little & H Bremer 1982
growthRate_RLB1982 = [0.4180,0.4470,0.4870,1.0790,1.1310,1.1810,1.2570,1.3130,1.2860,1.3590,1.2390,1.1890,1.8170,1.8640,1.9190]*log(2);
ppGpp_RLB1982      = [83.4290,73.2630,56.3370,40.7160,42.3690,35.2900,32.2710,31.8780,30.2380,30.5310,27.9370,25.9460,13.6520,11.6100,12.5670]/convFac;
plot(growthRate_RLB1982,ppGpp_RLB1982,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(2,:));

%   Paper 3: Kinetic properties of rrn promoters in Escherichia coli
%   X Zhang et al. (2002)
growthRate_Zhang2002 = [0.663,1.181,2.137,2.726]*log(2);
ppGpp_Zhang2002      = [51.535,33.029,12.926,3.909]/convFac;
plot(growthRate_Zhang2002,ppGpp_Zhang2002,'ko','MarkerSize',MS,'LineWidth',LW,'MarkerFaceColor',ColorPalette(3,:));

xlim([0,2]);
ylim([0,500]);

%----------------------------------------------------------------
%   Amino acid level vs Growth rate
%----------------------------------------------------------------

subplot(1,3,3);
hold on;

%   new marker size
MS = 10;

plot(growthRate,MET,'k-','LineWidth',LW);
axis square;
box on;
set(gca,'YScale','log');
xlabel('Growth rate (h^{-1})');
ylabel('Free amino acids (\muM)');

% %   Paper: Amino acid pool of Escherichia coli during the different phases of growth
% %   Growth rate calcated from growth curve in Fig. 4
% % growthrate_RR1970           = 0.447;    %   1/h
% % aminoacid_RR1970            = [168,231,41,41,6,10,9,45,7,638,171,3,8];  %   uM
% % plot(growthrate_RR1970,aminoacid_RR1970,'k.','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(1,:));
% % plot([growthrate_RR1970-0.05,growthrate_RR1970+0.05],[median(aminoacid_RR1970),median(aminoacid_RR1970)],'k-','Color',ColorPalette(1,:),'LineWidth',LW);
% plot(0.447,0.7,'ko','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(1,:));

% %   Paper: Discharging tRNAs: a tug of war between translation and detoxification in Escherichia coli
% growthrate_IAK2016_slow     = 0.96;     %   1/h
% aminoacid_IAK2016_slow      = [0.027,0.002,0.062,0.008,0.009,0.024,0.012,...
%     0.057,0.008,0.065,0.011,0.185]*1000; %   uM
% plot(growthrate_IAK2016_slow,aminoacid_IAK2016_slow,'k.','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(2,:));
% plot([growthrate_IAK2016_slow-0.05,growthrate_IAK2016_slow+0.05],[median(aminoacid_IAK2016_slow),median(aminoacid_IAK2016_slow)],'k-','Color',ColorPalette(2,:),'LineWidth',LW);
% 
% growthrate_IAK2016_fast     = 1.38;     %   1/h
% aminoacid_IAK2016_fast      = [1.21,0.32,0.21,0.36,0.13,0.29,0.25,0.27,...
%     4.00,0.34,0.77,0.21,3.00]*1000;  %   uM
% plot(growthrate_IAK2016_fast,aminoacid_IAK2016_fast,'k.','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(2,:));
% plot([growthrate_IAK2016_fast-0.05,growthrate_IAK2016_fast+0.05],[median(aminoacid_IAK2016_fast),median(aminoacid_IAK2016_fast)],'k-','Color',ColorPalette(2,:),'LineWidth',LW);

% %   Paper: Absolute metabolite concentrations and implied enzyme active site occupancy in Escherichia coli.
% growthrate_BDB2009_slow     = log(2)/(139/60);  %   1/h
% aminoacid_BDB2009_slow      = [8.79e+2,5.40e+2,7.35e+3,4.48e+4,3.06e+3,...
%     9.75e+1,3.41e+2/2,5.54e+2,6.59e+1,2.74e+1,...
%     3.61e+2,9.55e+1,2.05e+1,5.22e+1,1.07e+3];    %   uM
% plot(growthrate_BDB2009_slow,aminoacid_BDB2009_slow,'k.','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(3,:));
% plot([growthrate_BDB2009_slow-0.05,growthrate_BDB2009_slow+0.05],[median(aminoacid_BDB2009_slow),median(aminoacid_BDB2009_slow)],'k-','Color',ColorPalette(3,:),'LineWidth',LW);
errorbar(log(2)/(139/60),6.59e1,6.59e1-4.38e1,9.92e1-6.59e1,'ko','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');

% growthrate_BDB2009_interm   = log(2)/(89/60);   %   1/h
% aminoacid_BDB2009_interm    = [1.77e+3,9.70e+2,9.30e+3,1.49e+5,4.95e+3,...
%     1.75e+2,4.38e+2/2,7.62e+2,1.29e+2,4.21e+1,...
%     4.51e+2,1.50e+2,2.36e+1,8.74e+1,2.29e+3];    %   uM
% plot(growthrate_BDB2009_interm,aminoacid_BDB2009_interm,'k.','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(3,:));
% plot(growthrate_BDB2009_interm,Met,'ko','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(3,:));
% plot([growthrate_BDB2009_interm-0.05,growthrate_BDB2009_interm+0.05],[median(aminoacid_BDB2009_interm),median(aminoacid_BDB2009_interm)],'k-','Color',ColorPalette(3,:),'LineWidth',LW);
errorbar(log(2)/(89/60),1.29e2,1.29e2-9.87e1,1.68e2-1.29e2,'ko','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');

% growthrate_BDB2009_fast     = log(2)/(77/60);   %   1/h
% aminoacid_BDB2009_fast      = [2.55e+3,5.69e+2,5.11e+2,4.23e+3,9.60e+4,...
%     3.81e+3,6.76e+1,3.03e+2/2,4.05e+2,1.45e+2,...
%     1.82e+1,3.85e+2,6.80e+1,1.79e+2,1.21e+1,...
%     2.89e+1,4.02e+3];    %   uM
% plot(growthrate_BDB2009_fast,aminoacid_BDB2009_fast,'k.','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(3,:));
% plot(growthrate_BDB2009_fast,Met,'ko','MarkerSize',MS,'MarkerEdgeColor',ColorPalette(3,:));
% plot([growthrate_BDB2009_fast-0.05,growthrate_BDB2009_fast+0.05],[median(aminoacid_BDB2009_fast),median(aminoacid_BDB2009_fast)],'k-','Color',ColorPalette(3,:),'LineWidth',LW);
errorbar(log(2)/(77/60),1.45e2,1.45e2-1.31e2,1.61e2-1.45e2,'ko','MarkerSize',MS,'MarkerEdgeColor','k','MarkerFaceColor','w');

xlim([0,2]);
ylim([1e0,1e+5]);
clear all;

fnr_mutant = false;
arcA_mutant = false;

Time_Final = 25 * 3600;
GLC_ex_input = 10.0;

a1 = (0.01:0.02:0.98);
a2 = (1.0:1.0:15.0);
a3 = (20.0:10.0:50.0);
a = [a1 a2 a3];

Num_of_State_Variable = 36;
Num_of_Flux = 10;
options=odeset('RelTol',1e-5,'AbsTol',1e-20);
xx = zeros(length((0.0:60.0:Time_Final)),Num_of_State_Variable,length(a));
ff = zeros(length((0.0:60.0:Time_Final)),Num_of_Flux,length(a));

for i=1:length(a)
    
    IC = Initial_Concentration();
    called_by_ode15s = true;
    [t_0, x_0] = ode15s(@Kinetic_model,[-10.0*3600.0 0.0],IC,options,arcA_mutant,fnr_mutant,a(i),called_by_ode15s);
    [t, x] = ode15s(@Kinetic_model,(0.0:60.0:Time_Final),x_0(end,:),options,arcA_mutant,fnr_mutant,a(i),called_by_ode15s);

    called_by_ode15s = false;
    f = zeros(size(x,1),Num_of_Flux);
    for indx = 1:size(x,1)
        f(indx,:) = Kinetic_model(0,x(indx,:),arcA_mutant,fnr_mutant,a(i),called_by_ode15s);
    end
    
    if size(x,1) < size(xx(:,:,i),1)
        xx(:,:,i) = [x(1:end,:); zeros(length(xx(:,1,i))-length(x(:,1)),Num_of_State_Variable)];
        ff(:,:,i) = [f(1:end,:); zeros(length(xx(:,1,i))-length(f(:,1)),Num_of_Flux)];
    else
        xx(:,:,i) = x;
        ff(:,:,i) = f;
    end
    
end

t = (1.0/3600.0)*t;
ff = 3600.0*ff;

OD(:,1:length(a)) = xx(:,1,1:length(a));
GLC_ex(:,1:length(a)) = xx(:,2,1:length(a));
ACT_ex(:,1:length(a)) = xx(:,3,1:length(a));
LAC(:,1:length(a)) = xx(:,29,1:length(a));
ETH(:,1:length(a)) = xx(:,31,1:length(a));
FOR(:,1:length(a)) = xx(:,32,1:length(a));
Q(:,1:length(a)) = xx(:,34,1:length(a));
NADH(:,1:length(a)) = xx(:,35,1:length(a));
SUC_ex(:,1:length(a)) = xx(:,36,1:length(a));

TF_Fnr(:,1:length(a)) = ff(:,1,1:length(a))./3600.0;
TF_ArcA(:,1:length(a)) = ff(:,2,1:length(a))./3600.0;
v_Cyo(:,1:length(a)) = ff(:,3,1:length(a));
v_Cyd(:,1:length(a)) = ff(:,4,1:length(a));
v_PDH(:,1:length(a)) = ff(:,5,1:length(a));
v_Pfl(:,1:length(a)) = ff(:,6,1:length(a));
v_LDH(:,1:length(a)) = ff(:,7,1:length(a));
v_ADH(:,1:length(a)) = ff(:,8,1:length(a));
v_FRD(:,1:length(a)) = ff(:,9,1:length(a));
NAD(:,1:length(a)) = ff(:,10,1:length(a))./3600.0;

GLC_ex = real(GLC_ex);
Index = zeros(1,length(a));
for i=1:length(a)
    [IDX, D] = knnsearch(GLC_ex(:,i),0.01*GLC_ex(1,i));
    Index(1,i) = IDX;
end

GLC_zero_OD = zeros(1,length(a));
GLC_zero_GLC_ex = zeros(1,length(a));
GLC_zero_ACT_ex = zeros(1,length(a));
GLC_zero_LAC = zeros(1,length(a));
GLC_zero_ETH = zeros(1,length(a));
GLC_zero_FOR = zeros(1,length(a));
GLC_zero_Q = zeros(1,length(a));
GLC_zero_NADH = zeros(1,length(a));
GLC_zero_SUC_ex = zeros(1,length(a));
GLC_zero_TF_Fnr = zeros(1,length(a));
GLC_zero_TF_ArcA = zeros(1,length(a));
GLC_zero_v_Cyo = zeros(1,length(a));
GLC_zero_v_Cyd = zeros(1,length(a));
GLC_zero_v_PDH = zeros(1,length(a));
GLC_zero_v_Pfl = zeros(1,length(a));
GLC_zero_v_LDH = zeros(1,length(a));
GLC_zero_v_ADH = zeros(1,length(a));
GLC_zero_v_FRD = zeros(1,length(a));
GLC_zero_NAD = zeros(1,length(a));
for i=1:length(a)
    GLC_zero_OD(1,i) = OD(Index(i),i);
    GLC_zero_GLC_ex(1,i) = GLC_ex(Index(i),i);
    GLC_zero_ACT_ex(1,i) = ACT_ex(Index(i),i);
    GLC_zero_LAC(1,i) = LAC(Index(i),i);
    GLC_zero_ETH(1,i) = ETH(Index(i),i);
    GLC_zero_FOR(1,i) = FOR(Index(i),i);
    GLC_zero_Q(1,i) = Q(Index(i),i);
    GLC_zero_NADH(1,i) = NADH(Index(i),i);
    GLC_zero_SUC_ex(1,i) = SUC_ex(Index(i),i);

    GLC_zero_TF_Fnr(1,i) = TF_Fnr(Index(i),i);
    GLC_zero_TF_ArcA(1,i) = TF_ArcA(Index(i),i);
    GLC_zero_v_Cyo(1,i) = v_Cyo(Index(i),i);
    GLC_zero_v_Cyd(1,i) = v_Cyd(Index(i),i);
    GLC_zero_v_PDH(1,i) = v_PDH(Index(i),i);
    GLC_zero_v_Pfl(1,i) = v_Pfl(Index(i),i);
    GLC_zero_v_LDH(1,i) = v_LDH(Index(i),i);
    GLC_zero_v_ADH(1,i) = v_ADH(Index(i),i);
    GLC_zero_v_FRD(1,i) = v_FRD(Index(i),i);
    GLC_zero_NAD(1,i) = NAD(Index(i),i);
end

scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/8 scrsz(3)/5.5 scrsz(4)/4],'Name','Figure(1)')
plot(a, GLC_zero_ACT_ex,'r-','LineWidth',2.5); hold on
plot(a, GLC_zero_LAC,'b-','LineWidth',2.5); hold on
plot(a, GLC_zero_FOR,'c-','LineWidth',2.5); hold on
plot(a, GLC_zero_ETH,'g-','LineWidth',2.5); hold on
plot(a, GLC_zero_SUC_ex,'m-','LineWidth',2.5); hold on
legend('ACE','LAC','FOR','ETH','SUC');
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'xtick',(0:10:50),'ylim', [0 5.0],'ytick',(0:5));
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',12);
ylabel('Concentration [g/l]','fontname','arial','fontweight','bold','fontsize',12);
set(gcf,'color','white'); hold on

figure('Position',[1 scrsz(4)/8 scrsz(3)/2 scrsz(4)/1.5],'Name','Figure(2)')
b = [7.0 20.0];
subplot(3,3,1); plot(a, GLC_zero_TF_ArcA,'k-','LineWidth',2.5); hold on
subplot(3,3,1); plot(a, GLC_zero_TF_Fnr,'b-','LineWidth',2.5); hold on
subplot(3,3,1); plot([b(1) b(1) b(1)], [0.0 0.5 1.0],'k:','LineWidth',2.0); hold on
subplot(3,3,1); plot([b(2) b(2) b(2)], [0.0 0.5 1.0],'k:','LineWidth',2.0); hold on
subplot(3,3,1); plot(0,1,'ko','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8,'LineWidth',2.0); hold on
text(2.8,1.1,'II','FontSize',14)
text(11.8,1.1,'III','FontSize',14)
text(31.8,1.1,'IV','FontSize',14)
text(0.01,1.3,'I','FontSize',14)
text(0.01,1.15,'«','FontSize',14)
legend('TF_A_r_c_A','TF_F_n_r');
title('','fontname','arial','fontweight','bold','fontsize',12);
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 1.0],'ytick',[0 0.5 1.0]);
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('TF_F_n_r and TF_A_r_c_A [-]','fontname','arial','fontweight','bold','fontsize',10);

subplot(3,3,2); plot(a, GLC_zero_Q,'k-','LineWidth',2.5); hold on
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 2.0]);
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('Q [\mumol/gDCW]','fontname','arial','fontweight','bold','fontsize',10);
  
subplot(3,3,3); plot(a, GLC_zero_NADH./GLC_zero_NAD,'k-','LineWidth',2.5); hold on
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 2.0]);
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('NADH/NAD^+ [-]','fontname','arial','fontweight','bold','fontsize',10);
  
subplot(3,3,4); plot(a, (0.5*GLC_zero_v_Cyo+0.5*GLC_zero_v_Cyd)/1000.0,'k-','LineWidth',2.5); hold on
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 15.0]);
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('qOUR [mmol/gDCW/h]','fontname','arial','fontweight','bold','fontsize',10);

subplot(3,3,5); plot(a, GLC_zero_v_PDH/1000.0,'k-','LineWidth',2.5); hold on
subplot(3,3,5); plot(a, GLC_zero_v_Pfl/1000.0,'b-','LineWidth',2.5); hold on
legend('PDH','Pfl');
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 20.0]);
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('PDH and Pfl [mmol/gDCW/h]','fontname','arial','fontweight','bold','fontsize',10);
  
subplot(3,3,6); plot(a, GLC_zero_v_LDH/1000.0,'k-','LineWidth',2.5); hold on
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 5.0]);
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('LDH [mmol/gDCW/h]','fontname','arial','fontweight','bold','fontsize',10);
  
subplot(3,3,7); plot(a, GLC_zero_v_Cyo/1000.0,'k-','LineWidth',2.5); hold on
subplot(3,3,7); plot(a, GLC_zero_v_Cyd/1000.0,'b-','LineWidth',2.5); hold on
legend('Cyo','Cyd');
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 30.0]);
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('Cyo and Cyd [mmol/gDCW/h]','fontname','arial','fontweight','bold','fontsize',10);
   
subplot(3,3,8); plot(a, GLC_zero_v_ADH/1000.0,'k-','LineWidth',2.5); hold on
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 8.0],'ytick',(0:2:8));
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('ADH [mmol/gDCW/h]','fontname','arial','fontweight','bold','fontsize',10);
  
subplot(3,3,9); plot(a, GLC_zero_v_FRD/1000.0,'k-','LineWidth',2.5); hold on
set(gca,'fontname','arial','fontsize',10,'fontweight','bold','linewidth',1.5,'xlim',[0 a(end)],'ylim', [0 3.0],'ytick',(0:3));
xlabel('DO [%]','fontname','arial','fontweight','bold','fontsize',10);
ylabel('Frd [mmol/gDCW/h]','fontname','arial','fontweight','bold','fontsize',10);
set(gcf,'color','white'); hold on
  

 
    
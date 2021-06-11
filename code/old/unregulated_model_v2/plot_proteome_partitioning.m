clear all;
clc;

M=xlsread('msb145697-sup-0003-tables3.xlsx','Sheet2');

% C lim (1-5), A-lim (6-10), R-lim (11-14)
GrowthRate = log(2)./(reshape(M(:,1),14,6)/60); % C,A,R,U,S,O
GrowthRate = GrowthRate(:,1);
ProteomeFrac = reshape(M(:,2),14,6);

figure();

% plot ribosome vs nonribosome
NonR = sum(ProteomeFrac(:,[1,2,4,5,6]),2);
R = ProteomeFrac(:,3);

subplot(3,3,1);
hold on;
plot(GrowthRate(1:5),R(1:5),'ro-');
plot(GrowthRate(1:5),NonR(1:5),'b-o');
axis square;
box on;
axis([0.2,1.2,0,1]);

subplot(3,3,2);
hold on;
plot(GrowthRate(6:10),R(6:10),'ro-');
plot(GrowthRate(6:10),NonR(6:10),'b-o');
axis square;
box on;
axis([0.2,1.2,0,1]);

subplot(3,3,3);
hold on;
plot(GrowthRate(11:14),R(11:14),'ro-');
plot(GrowthRate(11:14),NonR(11:14),'b-o');
axis square;
box on;
axis([0.2,1.2,0,1]);

% plot ribosome, stress protein and the rest
R = ProteomeFrac(:,3);
S = ProteomeFrac(:,5);
E = sum(ProteomeFrac(:,[1,2,4,6]),2);

subplot(3,3,4);
hold on;
plot(GrowthRate(1:5),R(1:5),'ro-');
plot(GrowthRate(1:5),E(1:5),'b-o');
plot(GrowthRate(1:5),S(1:5),'k-o');
axis square;
box on;
axis([0.2,1.2,0,1]);

subplot(3,3,5);
hold on;
plot(GrowthRate(6:10),R(6:10),'ro-');
plot(GrowthRate(6:10),E(6:10),'b-o');
plot(GrowthRate(6:10),S(6:10),'k-o');
axis square;
box on;
axis([0.2,1.2,0,1]);

subplot(3,3,6);
hold on;
plot(GrowthRate(11:14),R(11:14),'ro-');
plot(GrowthRate(11:14),E(11:14),'b-o');
plot(GrowthRate(11:14),S(11:14),'k-o');
axis square;
box on;
axis([0.2,1.2,0,1]);

% plot R+S
RS = sum(ProteomeFrac(:,[3,5]),2);

subplot(3,3,7);
hold on;
plot(GrowthRate(1:5),RS(1:5),'ro-');
axis square;
box on;
axis([0.2,1.2,0,1]);

subplot(3,3,8);
hold on;
plot(GrowthRate(6:10),RS(6:10),'ro-');
axis square;
box on;
axis([0.2,1.2,0,1]);

subplot(3,3,9);
hold on;
plot(GrowthRate(11:14),RS(11:14),'ro-');
axis square;
box on;
axis([0.2,1.2,0,1]);

%% figure
AhpC = [1.4578	1.2129	1.0822	0.98587	1.0057, ...
        1.3454	1.2821	1.1905	0.95608	1,...
        0.79909	0.81241	0.88146	1];
KatG = [2.7343	1.8121	1.4379	1.1241	0.95036, ...
        1.5122	1.3508	1.2698	0.91347	1,...
        0.46675	0.50485	0.65448	1];

figure();

subplot(1,3,1);
hold on;
plot(GrowthRate(1:5),AhpC(1:5),'ro-');
plot(GrowthRate(1:5),KatG(1:5),'b-o');
axis square;
box on;
axis([0.2,1.2,0,3]);

subplot(1,3,2);
hold on;
plot(GrowthRate(6:10),AhpC(6:10),'ro-');
plot(GrowthRate(6:10),KatG(6:10),'b-o');
axis square;
box on;
axis([0.2,1.2,0,3]);

subplot(1,3,3);
hold on;
plot(GrowthRate(11:14),AhpC(11:14),'ro-');
plot(GrowthRate(11:14),KatG(11:14),'b-o');
axis square;
box on;
axis([0.2,1.2,0,3]);
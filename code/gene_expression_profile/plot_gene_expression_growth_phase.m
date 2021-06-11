clear all;
clc;

%% read table
gene_exp = readtable('./SuppTables_v2.xlsx');
gene_exp_array = table2array(gene_exp(:,2:end));

%% plot

figure();
time = [[0.6,4,20,40,60,90,120]/60,6.3,10.5,15,48]; % h first 0.6 should be zero

% Ribosomal genes
subplot(2,2,1);
hold on;

gnames = {'rpsT';'rplQ';'rpsD';'rpsK';'rpsM';'rpmJ';'secY';'rplO';'rpmD';'rpsE';'rplR';...
          'rplF';'rpsH';'rpsN';'rplE';'rplX';'rplN';'rpsQ';'rpmC';'rplP';'rpsC';'rplV';...
          'rpsS';'rplB';'rplW';'rplD';'rplC';'rpsJ';'tufA';'fusA';'rpsG';'rpsL';'rpoA';...
          'rpoB';'rpoC';'rpoZ'};
cls = jet(size(gnames,1));
for i=1:size(gnames,1)
    expression = gene_exp_array(find(strcmp(gene_exp.Var1,gnames{i})),:);
    plot(time,expression,'-','Color',cls(i,:));
end
legend(gnames);

set(gca,'Xscale','log');
set(gca,'Yscale','log');
axis([0.01,48,0.005,500]);

plot([4,4]/60,[0.01,100],'k-');
plot([120,120]/60,[0.01,100],'k-');
plot([4,120]/60,[100,100],'k-');
plot([4,120]/60,[0.01,0.01],'k-');

axis square;
box on;

xlabel('Time (h)');
ylabel('Relative gene expression');

% Metabolic proteins
subplot(2,2,2);
hold on;

gnames = {'acnA';'crp';'spoT'};
cls = jet(size(gnames,1));
for i=1:size(gnames,1)
    expression = gene_exp_array(find(strcmp(gene_exp.Var1,gnames{i})),:);
    plot(time,expression,'-','Color',cls(i,:));
end
legend(gnames);

set(gca,'Xscale','log');
set(gca,'Yscale','log');
axis([0.01,48,0.005,500]);

plot([4,4]/60,[0.01,100],'k-');
plot([120,120]/60,[0.01,100],'k-');
plot([4,120]/60,[100,100],'k-');
plot([4,120]/60,[0.01,0.01],'k-');

axis square;
box on;

xlabel('Time (h)');
ylabel('Relative gene expression');

% Stress proteins
subplot(2,2,3);
hold on;

gnames = {'dnaK';'dnaJ';'lon';'ahpC';'katG';'katE';'sodA';'sodB';'sodC';'rpoS'};
cls = jet(size(gnames,1));
for i=1:size(gnames,1)
    expression = gene_exp_array(find(strcmp(gene_exp.Var1,gnames{i})),:);
    plot(time,expression,'-','Color',cls(i,:));
end
legend(gnames);

set(gca,'Xscale','log');
set(gca,'Yscale','log');
axis([0.01,48,0.005,500]);

plot([4,4]/60,[0.01,100],'k-');
plot([120,120]/60,[0.01,100],'k-');
plot([4,120]/60,[100,100],'k-');
plot([4,120]/60,[0.01,0.01],'k-');

axis square;
box on;

xlabel('Time (h)');
ylabel('Relative gene expression');

% Division proteins
subplot(2,2,4);
hold on;

gnames = {'ftsZ';'ftsA';'ftsW';'ftsQ';'ftsX';'ftsE';'ftsY';'ftsN';'ddlB';...
          'cedA';'minD';'minE'};
cls = jet(size(gnames,1));
for i=1:size(gnames,1)
    expression = gene_exp_array(find(strcmp(gene_exp.Var1,gnames{i})),:);
    plot(time,expression,'-','Color',cls(i,:));
end
legend(gnames);

set(gca,'Xscale','log');
set(gca,'Yscale','log');
axis([0.01,48,0.005,500]);

plot([4,4]/60,[0.01,100],'k-');
plot([120,120]/60,[0.01,100],'k-');
plot([4,120]/60,[100,100],'k-');
plot([4,120]/60,[0.01,0.01],'k-');

axis square;
box on;

xlabel('Time (h)');
ylabel('Relative gene expression');
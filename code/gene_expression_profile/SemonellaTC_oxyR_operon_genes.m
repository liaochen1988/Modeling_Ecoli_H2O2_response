clear all;
clc;

%% read table
gene_exp = readtable('./SuppTables_v2.xlsx');
gene_exp_array = table2array(gene_exp(:,2:end));

%% plot

figure();
hold on;
time = [[0.6,4,20,40,60,90,120]/60,6.3,10.5,15,48]; % h first 0.6 should be zero

gnames = {'hemH';'ahpC';'fur';'dps';'grxA';...
          'ybjC';'hcp';'hcr';'poxB';'ychF';...
          'sufA';'sufB';'sufC';'sufD';'sufS';...
          'znuC';'znuA';'mntH';'gor';'metR';...
          'metE';'katG';'fhuF'};
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
clear all;
clc;

%% Fur regulon
% the first three values are protein synthesis rate in
% MOPS complete, MOPS minimal and MOPS complete w/o methionine
% the fourth value is aa. length
% the fifth value is molecular weight (kD)
% the sixth value is whether this gene is positively (1) or negatively(0)
% regulated

gene_names = {'Ipd';'fhuA';'fhuC';'fhuD';'fhuB';...
              'glnD';'cyoA';'cyoB';'cyoC';'cyoD';...
              'cyoE';'glnK';'amtB';'fepA';'entD';...
              'fes';'ybdZ';'entF';'fepE';'fepD';...
              'fepG';'fepC';'entS';'RutR';'fepB';...
              'entC';'entE';'entB';'entA';'entH';...
              'sdhC';'sdhD';'sdhA';'sdhB';'gpmA';...
              'fiu';'aspC';'ompF';'pyrC';'fhuE';...
              'ndh';'oppA';'oppB';'oppC';'oppD';...
              'oppF';'tonB';'fnr';'nohA';'ydfN';...
              'tfaQ';'fumC';'purR';'sufA';'sufB';...
              'sufC';'sufD';'sufS';'sufE';'katE';...
              'gdhA';'mntP';'zwf';'Fur';'flhD';...
              'flhC';'ftnA';'zinT';'nac';'gnd';...
              'rcnA';'rcnB';'cirA';'mntH';'hmp';...
              'grcA';'nrdH';'nrdI';'nrdE';'nrdF';...
              'rpoS';'exbD';'exbB';'nfeF';'garP';...
              'garL';'garR';'garK';'rnpB';'gltB';...
              'gltD';'gltF';'gspC';'gspD';'gspE';...
              'gspF';'gspG';'gspH';'gspI';'gspJ';...
              'gspK';'gspL';'gspM';'gspO';'feoA';...
              'feoB';'feoC';'ryhB';'yhhY';'sodA';...
              'metJ';'katG';'metH';'soxS';'soxR';...
              'fumB';'fecA';'fecB';'fecC';'fecD';...
              'fecE';'fecI';'fecR';'fhuF';'yjjZ'};
          
synthesis_rate = [476, 124, 254, 320, 35.844, 1;...              % hemH
                  92449, 37153, 76357, 187, 20.761, 1;...        % ahpC
                  5698,1554,4883,521,56.177,1;...              % ahpF
                  543,229,278,248,27.495,1;...                 % dsbG
                  275,119,132,28,3.109,1;...                  % uof
                  6426,2619,4625,148,16.795,1;...              % fur
                  2131,13008,2446,167,18.695,1;...             % dps
                  1247,431,807,85,9.685,1;...                 % grxA
                  46,38,28,95,10.521,0;...                     % ybjC
                  1916,646,1110,240,26.801,0;...               % nfsA
                  101,46,59,300,32.436,0;...                   % rimK
                  1365,367,745,158,17.666,0;...                % ybjN
                  4,1,0,550,60.064,1;...                       % hcp
                  18,16,0,322,35.74,1;...                     % hcr
                  40,877,24,572,62.011,1;...                   % poxB
                  11634,1907,5028,363,39.667,0;...             % ychF
                  146,359,126,122,13.3,1;...                 % sufA
                  75,143,71,495,54.745,1;...                   % sufB
                  117,226,94,248,27.582,1;...                  % sufC
                  82,152,60,423,46.823,1;...                   % sufD
                  76,117,23,138,15.8,1;...                   % sufE
                  63,122,33,406,44.434,1;...                   % sufS
                  332,468,1322,252,27.867,1;...                % znuC
                  252,301,450,261,27.729,1;...                 % znuB
                  3453,25238,57733,310,33.777,1;...            % znuA
                  360,20357,50070,216,24.762,1;...             % zinT
                  122,123,165,1039,106.82,0;...                % flu
                  136,58,57,412,44.194,1;...                   % mntH
                  1704,569,648,139,15.555,1;...                % trxC
                  3255,1055,1664,450,48.772,1;...              % gor
                  292,398,8396,317,35.629,1;...                % metR
                  571,63552,282818,753,84.673,1;...            % metE
                  2814,1908,1550,726,80.024,1;...              % katG
                  1867,726,1332,305,34.276,0;...               % oxyR
                  60,5,5,447,47.138,0;...                      % gntP
                  595,144,325,394,44.838,0;...                 % uxuA
                  364,90,213,486,53.58,0;...                  % uxuB
                  319,261,138,262,30.113,0;...                 % fhuF
                  25,21,2,78,8.697,1 ...                      % yjjZ
                  ];
                  
total_protein = [680,238,450]; % fg
doubling_time = [21.5,56.3,26.5]; % min

synthesis_rate_pos = synthesis_rate(find(synthesis_rate(:,end)==1),:);
synthesis_rate_neg = synthesis_rate(find(synthesis_rate(:,end)==0),:);
assert(size(synthesis_rate_pos,1)+size(synthesis_rate_neg,1)==size(synthesis_rate,1));

frac_pos = synthesis_rate_pos(:,1:3);
for i=1:size(frac_pos,1)
    % synthesis rate : molecules per generation
    % molecular weigth: in unit of kD
    % 1fg = 1e-15 g
    % 1 au = 1.66054e-24 g
    frac_pos(i,:) = frac_pos(i,:)*synthesis_rate_pos(i,5)*1e3*1.66054e-24*1e15./total_protein; % fg
end

frac_neg = synthesis_rate_neg(:,1:3);
for i=1:size(frac_neg,1)
    frac_neg(i,:) = frac_neg(i,:)*synthesis_rate_neg(i,5)*1e3*1.66054e-24*1e15./total_protein; % fg
end

gene_names_pos = gene_names(find(synthesis_rate(:,end)==1));
gene_names_neg = gene_names(find(synthesis_rate(:,end)==0));

%figure();

subplot(2,3,1);
hold on;
bar([1:3],frac_pos','stacked');
set(gca,'XTick',[1,2,3]);
set(gca,'XTicklabel',{'MOPS c';'MOPS m';'MOPS c w/o met'});
axis square;
box on;
legend(gene_names_pos);

subplot(2,3,4);
hold on;
bar([1:3],frac_neg','stacked');
set(gca,'XTick',[1,2,3]);
set(gca,'XTicklabel',{'MOPS c';'MOPS m';'MOPS c w/o met'});
axis square;
legend(gene_names_neg);

box on;
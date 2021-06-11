function hostparams = readParameters()

%   # of amino acids in R-sector and affiliated proteins (aa)
%   biological constant
hostparams.('Lr')            = 7336*1.6;

%   # of amino acids in MetE (aa)
hostparams.('Le')            = 753;

%   maximum peptide elongation rate (1/h)
%   close to 20 aa/s, the observed maximum elongation rate
hostparams.('kelong')        = 21*3600;

%   total amino acid concentration (uM)
%   Paper: Growth rate of Escherichia coli by A Marr (1991)
hostparams.('Pt')             = 3.00E6;

%   ppGpp decay rate (1/h)
%   Paper: Interaction of alleles of the relA, relC and spoT genes in
%   Escherichia coli: analysis of the interconversion of GTP, ppGpp, and
%   pppGpp by N Fiil et al. (1977)
hostparams.('kdPpGpp')          = 0.0231*3600;

%   maximum proteome allocation for R sector
%   Close to observed 0.55 in paper
%   Interdependence of cell growth and gene expression: origins and
%   consequences by M Scott et al. (2010)
hostparams.('phiQ')             = 0.52;

%   fraction of Methionine
%   biological constant
hostparams.('fmet')             = 0.027;

%   fraction of active ribosomes
%   Paper: Modulation of chemical composition and other parameters of the cell by growth rate
hostparams.('f_act')            = 0.80;

%   maximum ppGpp synthesis rate (1/h)
%   Paper: ?How fast-growing bacteria robustly tune their ribosome concentration to approximate growth-rate maximization
%   RelA concentration : 0.21 uM
hostparams.('ksPpGpp')          = 75*3600*0.21;

%-------------------------------
%   Free parameters
%-------------------------------

%   Michaelis constant of amino acid (uM)
hostparams.('Km_met')           = 40.00;

%   ppGpp threshold for ribosome synthesis (uM)
hostparams.('Ki_ppGpp')         = 60.00;

%   fraction of metE in all metablic proteins
%   Paper: Quantifying absolute protein synthesis rates reveals principles underlying allocation of cellular resources
%   This parameter is set to give growth rate = 1.47 at k_cat = 0.12
hostparams.('alpha_metE')       = 0.39;

end

